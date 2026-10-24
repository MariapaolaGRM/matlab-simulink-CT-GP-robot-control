%% Esegue controlli di consistenza sul dataset utilizzato per l'addestramento e la validazione del Gaussian Process
% X = [ddq1 ddq2 dq1 dq2 q1 q2],  Y = [tau_tilde1 tau_tilde2]

%clear; close all; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));
%% Caricamento
load(fullfile(root, 'data', 'dataset.mat'));

X  = D_train.X;   
Y  = D_train.Y;   

id = D_train.id;

Xv = D_val.X;     
Yv = D_val.Y;
ids = unique(id).';

%% Dimensioni, NaN/Inf, traiettorie
ok_dim = size(X,2)==6 && size(Y,2)==2 && size(X,1)==size(Y,1) ...
      && size(Xv,2)==6 && size(Yv,2)==2;
ok_nan = all(isfinite([X(:); Y(:); Xv(:); Yv(:)]));

fprintf('Training: %d punti, %d traiettorie\n', size(X,1), numel(ids));
fprintf('Validazione: %d punti, %d traiettorie\n', size(Xv,1), numel(unique(D_val.id)));

tol = 0.05;
ok_lim = all(X(:,5) >= q_lim(1,1)-tol & X(:,5) <= q_lim(1,2)+tol & ...
             X(:,6) >= q_lim(2,1)-tol & X(:,6) <= q_lim(2,2)+tol);
%% Coerenza analitica di tau_tilde
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

%% Coerenza tra q, dq, ddq (derivate numeriche per traiettoria)
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

%% Copertura delle traiettorie di test

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

% Traiettoria di test usata nel confronto (cubica 3 s)
[tte, Xte] = test_traj([pi/4; pi/8], [pi/2; pi/4], 3, 0.5, Ts_data);
d_te = nearest_distance((Xte - mu)./sg, Z);
pct  = 100*mean(d_te > soglia);
ok_cov = pct < 25;
fprintf('Copertura cubica 3 s: mediana %.2f | max %.2f | oltre soglia %.1f %%\n', ...
    median(d_te), max(d_te), pct);

d_val = nearest_distance((Xv - mu)./sg, Z);
fprintf('%-26s mediana %.2f | max %.2f | oltre soglia %.1f %%\n', ...
    'validazione', median(d_val), max(d_val), 100*mean(d_val > soglia));

%% Esito
fprintf('\n--- ESITO ---\n');
fprintf('Dimensioni X/Y:                    %s\n', esito(ok_dim));
fprintf('Nessun NaN/Inf:                    %s\n', esito(ok_nan));
fprintf('Configurazioni dentro q_lim:       %s\n', esito(ok_lim));
fprintf('Coerenza analitica di tau_tilde:   %s\n', esito(ok_ana));
fprintf('Coerenza q, dq, ddq:               %s\n', esito(ok_kin));
fprintf('Copertura traiettoria di test:     %s\n', esito(ok_cov));

if ok_dim && ok_nan && ok_lim && ok_ana && ok_kin && ok_cov
    fprintf('\nDataset pronto per il training del GP.\n');
else
    fprintf('\nCi sono controlli da risolvere prima di addestrare il GP.\n');
end

%% FUNZIONI
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

function [t, Xte] = test_traj(q0, qf, T, t_tail, Ts)
    % ingresso del GP lungo il riferimento cubico: [ddq_d dq_d q_d], con coda a riposo
    t    = (0:Ts:T+t_tail)';
    tt   = min(t,T);
    rest = double(t < T);
    Xte  = zeros(numel(t),6);
    for j = 1:2
        D  = qf(j) - q0(j);
        a2 = 3*D/T^2;   a3 = -2*D/T^3;
        Xte(:,j)   = (2*a2 + 6*a3*tt) .* rest;                 % accelerazione
        Xte(:,j+2) = (2*a2*tt + 3*a3*tt.^2) .* rest;           % velocità
        Xte(:,j+4) = q0(j) + a2*tt.^2 + a3*tt.^3;              % posizione
    end
end
