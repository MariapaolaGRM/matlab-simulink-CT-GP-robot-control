% Valuta le prestazioni dei Gaussian Process addestrati utilizzando
% un dataset di validazione indipendente da quello di addestramento.
clear; close all; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));

load('dataset_gp.mat');            % D_train, D_val
load('gp_trained.mat');            % gp, Xm, Ym (GP con m = 400)

Xv = D_val.X;  Yv = D_val.Y;  tv = D_val.t;  idv = D_val.id;

%% ===== PUNTO 1: errore per singola traiettoria di validazione =====
Yp = zeros(size(Yv));  Ysd = zeros(size(Yv));
for i = 1:2
    [Yp(:,i), Ysd(:,i)] = predict(gp{i}, Xv);
end

ids = unique(idv).';
nT  = numel(ids);

% distanza dai punti di training effettivamente usati dal GP (spazio normalizzato)
mu = mean(Xm);  sg = std(Xm);  sg(sg==0) = 1;
dmin = min(pdist2((Xv-mu)./sg, (Xm-mu)./sg), [], 2);

% tabella riassuntiva
rows = zeros(nT, 9);
for k = 1:nT
    in = idv == ids(k);
    e  = Yp(in,:) - Yv(in,:);                 % errore di predizione
    rows(k,:) = [ids(k), ...
        sqrt(mean(Yv(in,:).^2)), ...          % RMSE senza GP, predizione nulla(giunto 1, 2)
        sqrt(mean(e.^2)), ...                 % RMSE con GP   (giunto 1, 2)
        max(abs(Xv(in,3:4)),[],1), ...        % |dq| massima  (giunto 1, 2)
        max(abs(Xv(in,1:2)),[],1)];           % |ddq| massima (giunto 1, 2)
end
% la colonna 2-3 è senza GP, 4-5 con GP, 6-7 dq, 8-9 ddq
Tab = array2table(rows, 'VariableNames', ...
    {'traj','RMSE0_1','RMSE0_2','RMSE_1','RMSE_2','dqmax1','dqmax2','ddqmax1','ddqmax2'});
Tab.dist_max = arrayfun(@(k) max(dmin(idv==ids(k))), ids).';
Tab.dist_max = Tab.dist_max(:);
disp(Tab);

% grafici: una figura per giunto, un pannello per traiettoria, asse in secondi
for i = 1:2
    figure('Name', sprintf('Giunto %d', i));
    for k = 1:nT
        in = idv == ids(k);
        subplot(2, ceil(nT/2), k);
        plot(tv(in), Yv(in,i), 'k'); hold on
        plot(tv(in), Yp(in,i), 'r');
        plot(tv(in), Yp(in,i)+2*Ysd(in,i), 'r:');
        plot(tv(in), Yp(in,i)-2*Ysd(in,i), 'r:');
        grid on; xlabel('t (s)'); ylabel(sprintf('\\tau~_%d (N m)', i));
        title(sprintf('traj %d | RMSE %.1f (senza GP %.1f)', ids(k), ...
              Tab.(sprintf('RMSE_%d',i))(k), Tab.(sprintf('RMSE0_%d',i))(k)));
    end
    legend('vero','GP','\pm2\sigma','Location','best');
end

% relazione tra errore e velocità/distanza dal training
figure;
subplot(1,2,1);
scatter(Tab.ddqmax1, Tab.RMSE_1, 60, 'filled'); hold on
scatter(Tab.ddqmax2, Tab.RMSE_2, 60, 'filled'); grid on
xlabel('|ddq| massima (rad/s^2)'); ylabel('RMSE con GP (N m)'); legend('giunto 1','giunto 2');
subplot(1,2,2);
scatter(Tab.dist_max, Tab.RMSE_1, 60, 'filled'); hold on
scatter(Tab.dist_max, Tab.RMSE_2, 60, 'filled'); grid on
xlabel('distanza massima dal training (norm.)'); ylabel('RMSE con GP (N m)');

%% ===== PUNTO 2: curva di apprendimento =====
X = D_train.X;  Y = D_train.Y;
m_list = [25 50 100 200 400 800];
n_rep  = 5; % numero di volte in cui ripeto l’addestramento per ogni sottoinsieme di campioni
rmse_lc = nan(numel(m_list), n_rep, 2);
rmse0   = sqrt(mean(Yv.^2));                  % riferimento: predire sempre 0

rng(10);
for a = 1:numel(m_list)
    m = min(m_list(a), size(X,1));
    for r = 1:n_rep
        idx = randperm(size(X,1), m);
        for i = 1:2
            g = fitrgp(X(idx,:), Y(idx,i), ...
                'BasisFunction','none', ...
                'KernelFunction','squaredexponential', ...
                'Standardize',true, ...
                'FitMethod','exact','PredictMethod','exact');
            yp = predict(g, Xv);
            rmse_lc(a,r,i) = sqrt(mean((yp - Yv(:,i)).^2));
        end
    end
    fprintf('m = %4d | RMSE medio: [%.2f %.2f] N m\n', m, ...
        mean(rmse_lc(a,:,1)), mean(rmse_lc(a,:,2)));
end

figure;
for i = 1:2
    subplot(1,2,i);
    mu_lc = mean(rmse_lc(:,:,i), 2);
    sd_lc = std(rmse_lc(:,:,i), 0, 2);
    errorbar(m_list, mu_lc, sd_lc, 'o-', 'LineWidth', 1.5); hold on
    yline(rmse0(i), 'k--', 'senza GP');
    set(gca,'XScale','log'); grid on
    xlabel('numero di punti di training m'); ylabel('RMSE validazione (N m)');
    title(sprintf('Giunto %d', i));
end

%% Curva di apprendimento del GP - grafico unico
figure;
hold on;
grid on;
box on;

% Colori dei due giunti
c1 = [0.00 0.45 0.74];   % blu
c2 = [0.90 0.35 0.10];   % arancione

% Statistiche del RMSE sui diversi esperimenti
mu1 = mean(rmse_lc(:,:,1), 2);
sd1 = std(rmse_lc(:,:,1), 0, 2);

mu2 = mean(rmse_lc(:,:,2), 2);
sd2 = std(rmse_lc(:,:,2), 0, 2);

% Curve con barre d'errore
errorbar(m_list, mu1, sd1, 'o-', ...
    'Color', c1, 'MarkerFaceColor', c1, ...
    'LineWidth', 1.5, 'MarkerSize', 7);

errorbar(m_list, mu2, sd2, 's-', ...
    'Color', c2, 'MarkerFaceColor', c2, ...
    'LineWidth', 1.5, 'MarkerSize', 7);

% Baseline senza GP
yline(rmse0(1), '--', ...
    sprintf('Senza GP, giunto 1 (%.1f)', rmse0(1)), ...
    'Color', c1, 'LineWidth', 1.5);

yline(rmse0(2), '--', ...
    sprintf('Senza GP, giunto 2 (%.1f)', rmse0(2)), ...
    'Color', c2, 'LineWidth', 1.5);

% Formattazione
set(gca, 'XScale', 'linear', 'FontSize', 12);

xlabel('Numero di punti di training m');
ylabel('RMSE di validazione (N m)');
title('Curva di apprendimento del Gaussian Process');

legend('Giunto 1', 'Giunto 2', ...
    'Location', 'northeast');

xlim([min(m_list)-10, max(m_list)+10]);

save('analisi_gp_risultati.mat','Tab','rmse_lc','m_list','rmse0');