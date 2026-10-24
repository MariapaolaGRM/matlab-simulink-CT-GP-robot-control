% Prepara i parametri dei Gaussian Process addestrati per il loro utilizzo
% all'interno del modello Simulink del controllore.
%clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));

% Caricamento
load(fullfile(root, 'data', 'dataset.mat'));
load(fullfile(root, 'data', 'gp_trained.mat'));

% Parametri dei disturbi
k_dist = 0;  k_noise = 0;          % interruttori
t_d    = 1.5;                      % istante del gradino di coppia (s)
d_amp  = [40; 15];                 % ampiezza del disturbo (N m), circa 5-10% dei picchi di coppia
sig_q  = 1e-3;  sig_dq = 1e-2;     % deviazioni standard del rumore (rad, rad/s)
Ts_n   = 1e-3;  seed   = 1;

% Standardizzazione degli ingressi
m = size(Xm,1);
P.mu = gp{1}.PredictorLocation(:).';       % 1x6 (standardizzazione degli ingressi)
P.sg = gp{1}.PredictorScale(:).';
P.Z  = (Xm - P.mu)./P.sg;                  % m x 6, punti di training standardizzati

% Estrazione parametri del GP
P.alpha = zeros(m,2);  
P.sL = zeros(1,2);  % la scala caratteristica del kernel, cioè quanto rapidamente cambia la correlazione tra punti
P.sF = zeros(1,2);  % l'ampiezza caratteristica della funzione appresa

D2 = pdist2(P.Z, P.Z).^2;
for i = 1:2
    kp = gp{i}.KernelInformation.KernelParameters;   % [lambda; sigma_f]
    P.sL(i) = kp(1);   P.sF(i) = kp(2);
    sn = gp{i}.Sigma; % il livello di rumore stimato dal modello
    K  = P.sF(i)^2 * exp(-0.5*D2/P.sL(i)^2);
    P.alpha(:,i) = (K + sn^2*eye(m)) \ Ym(:,i);      % alpha = (K + sn^2 I)^-1 y
end

% Verifica contro predict
Xv = D_val.X;
f1 = zeros(size(Xv,1),2);  
f2 = f1;
for n = 1:size(Xv,1), f1(n,:) = gp_mean(Xv(n,:).', P).'; end
for i = 1:2, f2(:,i) = predict(gp{i}, Xv); end
err = max(abs(f1 - f2),[],1);
fprintf('Scarto massimo gp_mean vs predict: [%.2e %.2e] N m\n', err);
if any(err > 1e-4), warning('Scarto troppo grande: controlla le convenzioni.'); end

save('gp_params.mat','-struct','P');       % ogni campo diventa una variabile