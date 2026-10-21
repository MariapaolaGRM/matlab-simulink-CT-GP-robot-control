% Esegue controlli di consistenza sul dataset utilizzato per l'addestramento
% e la validazione del Gaussian Process.
% X = [ddq1 ddq2 dq1 dq2 q1 q2],  Y = [tau_tilde1 tau_tilde2]
clear; close all; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));
%% 1) Caricamento
load('dataset_gp.mat');            % D_train, D_val, Ts_data, q_lim
X  = D_train.X;   
Y  = D_train.Y;   

id = D_train.id;

Xv = D_val.X;     
Yv = D_val.Y;
ids = unique(id).';

%% 2) Dimensioni, NaN/Inf, traiettorie
ok_dim = size(X,2)==6 && size(Y,2)==2 && size(X,1)==size(Y,1) ...
      && size(Xv,2)==6 && size(Yv,2)==2;
ok_nan = all(isfinite([X(:); Y(:); Xv(:); Yv(:)]));

ty   = D_train.types(ids);
npts = arrayfun(@(k) sum(id==k), ids);
fprintf('Training: %d punti, %d traiettorie (%d quintiche, %d trapezoidali, %d cubiche)\n', ...
    size(X,1), numel(ids), sum(strcmp(ty,'quintic')), ...
    sum(strcmp(ty,'trapezoid')), sum(strcmp(ty,'cubic')));
fprintf('Punti per traiettoria: min %d | mediana %d | max %d\n', ...
    min(npts), median(npts), max(npts));
fprintf('Validazione: %d punti, %d traiettorie\n', size(Xv,1), numel(unique(D_val.id)));

%% 3) Intervalli delle variabili
nomi = {'ddq1','ddq2','dq1','dq2','q1','q2','tau~1','tau~2'};
V = [X Y];
fprintf('\n%-8s %10s %10s %10s\n','','min','max','std');
for i = 1:8
    fprintf('%-8s %10.3f %10.3f %10.3f\n', nomi{i}, min(V(:,i)), max(V(:,i)), std(V(:,i)));
end
tol = 0.05;
ok_lim = all(X(:,5) >= q_lim(1,1)-tol & X(:,5) <= q_lim(1,2)+tol & ...
             X(:,6) >= q_lim(2,1)-tol & X(:,6) <= q_lim(2,2)+tol);

%% 4) Coerenza analitica di tau_tilde
% tau_tilde = tau_reale(ddq,dq,q) - tau_hat(ddq,dq,q), con gli stessi stati misurati
pr = get_params(true);   ph = get_params(false);
Yan = zeros(size(Y));
for n = 1:size(X,1)
    ddq = X(n,1:2).';   dq = X(n,3:4).';   q = X(n,5:6).';
    tau  = H(q,pr)*ddq + C(q,dq,pr)*dq + g(q,pr) + pr.Fv.*dq;
    tauh = H(q,ph)*ddq + C(q,dq,ph)*dq + g(q,ph);
    Yan(n,:) = (tau - tauh).';
end
err_ana = max(abs(Y - Yan),[],'all');
ok_ana  = err_ana < 1e-6;
fprintf('\nmax|Y - Y_analitico| = %.2e N m\n', err_ana);

%% 5) Coerenza tra q, dq, ddq (derivate numeriche per traiettoria)
err_kin = zeros(numel(ids),2);
for k = 1:numel(ids)
    in  = find(id == ids(k));
    q   = X(in,5:6);   dq = X(in,3:4);   ddq = X(in,1:2);
    dq_num  = (q(3:end,:)  - q(1:end-2,:))  / (2*Ts_data);
    ddq_num = (dq(3:end,:) - dq(1:end-2,:)) / (2*Ts_data);
    err_kin(k,1) = median(abs(dq_num  - dq(2:end-1,:)),'all')  / max(max(abs(dq),[],'all'),eps);
    err_kin(k,2) = median(abs(ddq_num - ddq(2:end-1,:)),'all') / max(max(abs(ddq),[],'all'),eps);
end
ok_kin = all(err_kin(:) < 0.1);
fprintf('Errore relativo mediano (dq da q | ddq da dq): max %.3f | %.3f\n', max(err_kin));

%% 6) Segni delle velocita e errore di inseguimento
mov = abs(X(:,3:4)) > 1e-3;
frac_pos = sum(X(:,3:4) > 1e-3) ./ sum(mov);
ok_sign  = all(frac_pos > 0.35 & frac_pos < 0.65);
fprintf('Quota di dq > 0 (tra i campioni in moto): giunto 1 %.0f %% | giunto 2 %.0f %%\n', ...
    100*frac_pos);

ok_trk = true;
if isfield(D_train,'e')
    e = D_train.e;
    fprintf('Errore di inseguimento in raccolta: RMS [%.4f %.4f] rad | max [%.4f %.4f] rad\n', ...
        sqrt(mean(e.^2)), max(abs(e)));
    ok_trk = max(abs(e),[],'all') < 0.05;
else
    fprintf('Errore di inseguimento: non disponibile (aggiungi D.e in collect_dataset)\n');
end

%% 7) Copertura delle traiettorie di test
tests = { 'cubica 3 s',               'cubic', [pi/4; pi/8], [pi/2; pi/4], 3;
          'trapezoidale corso 5 s',   'trap',  [0; pi/4],    [pi/2; pi/2], 5;
          'cubica veloce 1.5 s',      'cubic', [pi/4; pi/8], [pi/2; pi/4], 1.5 };
t_tail = 0.5;
nt  = size(tests,1);
col = lines(nt);

% normalizzazione con il solo training
mu = mean(X,1);   sg = std(X,0,1);   sg(sg==0) = 1;
Z  = (X - mu)./sg;

% riferimento: distanza di ogni traiettoria dalle ALTRE traiettorie
d_ref = zeros(size(Z,1),1);
for k = ids
    in = id == k;
    d_ref(in) = nearest_distance(Z(in,:), Z(~in,:));
end
s = sort(d_ref);   soglia = s(ceil(0.9*numel(s)));
fprintf('\nSoglia (90 perc. del riferimento): %.2f | mediana riferimento %.2f\n', soglia, median(d_ref));

Xte = cell(1,nt);   tte = cell(1,nt);   d_te = cell(1,nt);   ok_cov = false(1,nt);
for i = 1:nt
    [tte{i}, Xte{i}] = test_traj(tests{i,2}, tests{i,3}, tests{i,4}, tests{i,5}, t_tail, Ts_data);
    d_te{i} = nearest_distance((Xte{i} - mu)./sg, Z);
    pct = 100*mean(d_te{i} > soglia);
    [dmax, imax] = max(d_te{i});
    ok_cov(i) = pct < 25;
    fprintf('%-26s mediana %.2f | max %.2f (t = %.2f s) | oltre soglia %.1f %%\n', ...
        tests{i,1}, median(d_te{i}), dmax, tte{i}(imax), pct);
end
d_val = nearest_distance((Xv - mu)./sg, Z);
fprintf('%-26s mediana %.2f | max %.2f | oltre soglia %.1f %%\n', ...
    'validazione', median(d_val), max(d_val), 100*mean(d_val > soglia));

% grafici
figure('Name','Copertura');
subplot(1,3,1); plot(X(:,5),X(:,6),'.','Color',[.75 .75 .75]); hold on
for i = 1:nt, plot(Xte{i}(:,5),Xte{i}(:,6),'.-','Color',col(i,:)); end
grid on; xlabel('q_1 (rad)'); ylabel('q_2 (rad)'); title('configurazioni')
subplot(1,3,2); plot(X(:,3),X(:,1),'.','Color',[.75 .75 .75]); hold on
for i = 1:nt, plot(Xte{i}(:,3),Xte{i}(:,1),'.-','Color',col(i,:)); end
grid on; xlabel('dq_1 (rad/s)'); ylabel('ddq_1 (rad/s^2)'); title('giunto 1')
subplot(1,3,3); plot(X(:,4),X(:,2),'.','Color',[.75 .75 .75]); hold on
for i = 1:nt, plot(Xte{i}(:,4),Xte{i}(:,2),'.-','Color',col(i,:)); end
grid on; xlabel('dq_2 (rad/s)'); ylabel('ddq_2 (rad/s^2)'); title('giunto 2')
legend([{'training'}; tests(:,1)], 'Location','best');

figure('Name','Distanza dal training'); hold on
for i = 1:nt, plot(tte{i}, d_te{i}, 'Color', col(i,:), 'LineWidth', 1.5); end
yline(soglia,'k--','soglia');
grid on; xlabel('t (s)'); ylabel('distanza normalizzata dal training piu vicino');
legend([tests(:,1); {'soglia'}], 'Location','best');

%% 8) Esito
fprintf('\n--- ESITO ---\n');
fprintf('Dimensioni X/Y:                    %s\n', esito(ok_dim));
fprintf('Nessun NaN/Inf:                    %s\n', esito(ok_nan));
fprintf('Configurazioni dentro q_lim:       %s\n', esito(ok_lim));
fprintf('Coerenza analitica di tau_tilde:   %s\n', esito(ok_ana));
fprintf('Coerenza q, dq, ddq:               %s\n', esito(ok_kin));
fprintf('Segni di dq bilanciati:            %s\n', esito(ok_sign));
fprintf('Errore di inseguimento contenuto:  %s\n', esito(ok_trk));
for i = 1:nt
    fprintf('Copertura %-26s %s\n', [tests{i,1} ':'], esito(ok_cov(i)));
end
if ok_dim && ok_nan && ok_lim && ok_ana && ok_kin && ok_sign && ok_trk && all(ok_cov)
    fprintf('\nDataset pronto per il training del GP.\n');
else
    fprintf('\nCi sono controlli da risolvere prima di addestrare il GP.\n');
end

%% Funzioni locali
function s = esito(c)
    if c, s = 'OK'; else, s = 'ATTENZIONE'; end
end

function dmin = nearest_distance(A, B)
    % distanza euclidea di ogni riga di A dal punto piu vicino di B
    dmin = zeros(size(A,1),1);
    for i = 1:size(A,1)
        dmin(i) = sqrt(min(sum((B - A(i,:)).^2, 2)));
    end
end

function [t, Xte] = test_traj(type, q0, qf, T, t_tail, Ts)
    % ingresso del GP lungo il riferimento: [ddq_d dq_d q_d], con coda a riposo
    t    = (0:Ts:T+t_tail)';
    tt   = min(t,T);
    rest = double(t < T);
    Xte  = zeros(numel(t),6);
    for j = 1:2
        D = qf(j) - q0(j);
        switch type
            case 'cubic'
                a2 = 3*D/T^2;   a3 = -2*D/T^3;
                q = q0(j) + a2*tt.^2 + a3*tt.^3;
                v = 2*a2*tt + 3*a3*tt.^2;
                a = 2*a2 + 6*a3*tt;
            case 'trap'                       % rampe di durata T/3 (come nella relazione del corso)
                tc = T/3;   vm = D/(T - tc);   ac = vm/tc;
                i1 = tt < tc;   i3 = tt > T - tc;   i2 = ~i1 & ~i3;
                q = zeros(size(t));   v = q;   a = q;
                q(i1) = q0(j) + 0.5*ac*tt(i1).^2;
                q(i2) = q0(j) + vm*(tt(i2) - tc/2);
                q(i3) = qf(j) - 0.5*ac*(T - tt(i3)).^2;
                v(i1) = ac*tt(i1);   v(i2) = vm;   v(i3) = ac*(T - tt(i3));
                a(i1) = ac;          a(i3) = -ac;
            otherwise
                error('Tipo di test sconosciuto: %s', type);
        end
        Xte(:,j)   = a .* rest;
        Xte(:,j+2) = v .* rest;
        Xte(:,j+4) = q;
    end
end