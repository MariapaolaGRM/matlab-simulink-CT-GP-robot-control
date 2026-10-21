% Addestra due modelli Gaussian Process indipendenti, uno per ciascun giunto
% del manipolatore robotico, utilizzando il dataset generato in precedenza.
clear; close all; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));

load('dataset_gp.mat');                 % D_train, D_val
root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));
X  = D_train.X;                         % [ddq1 ddq2 dq1 dq2 q1 q2]
Y  = D_train.Y;                         % [tau_tilde1 tau_tilde2]
Xv = D_val.X;   Yv = D_val.Y;

%% 1) Sottocampionamento a m punti
m = 400;
rng(3);
idx = randperm(size(X,1), min(m,size(X,1)));
Xm = X(idx,:);   Ym = Y(idx,:);

%% 2) Un GP per giunto (media zero, kernel SE, verosimiglianza)
gp = cell(1,2);
for i = 1:2
    gp{i} = fitrgp(Xm, Ym(:,i), ...
        'BasisFunction','none', ...
        'KernelFunction','squaredexponential', ...
        'Standardize',true, ...
        'FitMethod','exact','PredictMethod','exact');
    fprintf('Giunto %d: sigma_n = %.3g | lambda = %.3g | sigma_f = %.3g\n', i, ...
        gp{i}.Sigma, gp{i}.KernelInformation.KernelParameters(1), ...
        gp{i}.KernelInformation.KernelParameters(2));
end

%% 3) Validazione sulle cubiche
Yp = zeros(size(Yv));  Ysd = zeros(size(Yv));
for i = 1:2
    [Yp(:,i), Ysd(:,i)] = predict(gp{i}, Xv);
end
rmse   = sqrt(mean((Yp - Yv).^2));
rmse0  = sqrt(mean(Yv.^2));             % errore se il GP predicesse sempre 0
fprintf('RMSE validazione: [%.3f %.3f] N*m | senza GP: [%.3f %.3f] N*m\n', rmse, rmse0);

figure;
for i = 1:2
    subplot(2,1,i);
    plot(Yv(:,i),'k'); hold on; plot(Yp(:,i),'r');
    plot(Yp(:,i)+2*Ysd(:,i),'r:'); plot(Yp(:,i)-2*Ysd(:,i),'r:');
    grid on; ylabel(sprintf('\\tau~_%d (N m)',i));
    legend('vero','GP','\pm2\sigma');
end

save('gp_trained.mat','gp','Xm','Ym');