% Confronta le prestazioni di diverse strategie di controllo del manipolatore
% robotico a due giunti mediante simulazioni in MATLAB/Simulink.
clear; close all; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));

model = 'CT_PD_robot_2R';
Ts = 1e-3;  tf = 3;  Tsim = 6;

%% 1) Guadagni LG e HG
p_hat = get_params(false);
q2s = linspace(-pi,pi,200);  Hd = zeros(2,numel(q2s));
for k = 1:numel(q2s), Hd(:,k) = diag(H([0;q2s(k)],p_hat)); end
H_nom = min(Hd,[],2);  zeta = 1;
Kp_LG = diag(H_nom*4^2);   Kd_LG = diag(2*zeta*4*H_nom);
Kp_HG = diag(H_nom*20^2);  Kd_HG = diag(2*zeta*20*H_nom);

%% 2) Traiettoria di test (cubica, tenuta ferma dopo tf)
q0 = [pi/4; pi/8];  qf = [pi/2; pi/4];  dq0 = [0;0];
t  = (0:Ts:Tsim)';  tc = min(t,tf);  on = double(t <= tf);
q_des = zeros(numel(t),2);  qd_des = q_des;  qdd_des = q_des;
for j = 1:2
    D = qf(j)-q0(j);  a2 = 3*D/tf^2;  a3 = -2*D/tf^3;
    q_des(:,j)   = q0(j) + a2*tc.^2 + a3*tc.^3;
    qd_des(:,j)  = (2*a2*t + 3*a3*t.^2).*on;
    qdd_des(:,j) = (2*a2 + 6*a3*t).*on;
end

%% 2) Traiettoria di test (trapezoidale, durata 5 s)
% q0 = [pi/4; pi/8];
% qf = [pi/2; pi/4];
% dq0 = [0; 0];
% 
% %tf = 5;           % durata del movimento [s]
% Tsim = 6;          % un secondo di mantenimento finale
% t = (0:Ts:Tsim)';
% 
% q_des   = zeros(numel(t),2);
% qd_des  = zeros(numel(t),2);
% qdd_des = zeros(numel(t),2);
% 
% % Profilo trapezoidale: 30% accelerazione, 40% velocità
% % costante, 30% decelerazione
% t_acc = 0.30*tf;
% t_dec = tf - t_acc;
% tc_frac = 0.30;
% 
% for j = 1:2
%     D = qf(j) - q0(j);
% 
%     % Accelerazione costante e velocità massima
%     acc = D/(t_acc*(tf-t_acc));
%     v_max = acc*t_acc;
% 
%     for k = 1:numel(t)
%         tk = t(k);
% 
%         if tk <= t_acc
%             % Fase 1: accelerazione
%             q_des(k,j) = q0(j) + 0.5*acc*tk^2;
%             qd_des(k,j) = acc*tk;
%             qdd_des(k,j) = acc;
% 
%         elseif tk <= t_dec
%             % Fase 2: velocità costante
%             q_des(k,j) = q0(j) + ...
%                 0.5*acc*t_acc^2 + v_max*(tk-t_acc);
%             qd_des(k,j) = v_max;
%             qdd_des(k,j) = 0;
% 
%         elseif tk <= tf
%             % Fase 3: decelerazione
%             tr = tf-tk;
%             q_des(k,j) = qf(j) - 0.5*acc*tr^2;
%             qd_des(k,j) = acc*tr;
%             qdd_des(k,j) = -acc;
% 
%         else
%             % Dopo 5 s: posizione finale mantenuta
%             q_des(k,j) = qf(j);
%             qd_des(k,j) = 0;
%             qdd_des(k,j) = 0;
%         end
%     end
% end

% Visualizzazione della traiettoria desiderata
figure;
for j = 1:2
    subplot(3,1,1);
    plot(t,q_des(:,j)); hold on; grid on;
    ylabel('q_d (rad)');

    subplot(3,1,2);
    plot(t,qd_des(:,j)); hold on; grid on;
    ylabel('dq_d (rad/s)');

    subplot(3,1,3);
    plot(t,qdd_des(:,j)); hold on; grid on;
    ylabel('ddq_d (rad/s^2)');
    xlabel('t (s)');
end
%% 3) Quattro controllori
names = {'LG-PD','HG-PD','CT','CT-GP'};
cfg = [struct('Kp',Kp_LG,'Kd',Kd_LG,'k_ct',0,'k_gp',0), ...
       struct('Kp',Kp_HG,'Kd',Kd_HG,'k_ct',0,'k_gp',0), ...
       struct('Kp',Kp_LG,'Kd',Kd_LG,'k_ct',1,'k_gp',0), ...
       struct('Kp',Kp_LG,'Kd',Kd_LG,'k_ct',1,'k_gp',1)];

load_system(model);
set_param(model,'SolverType','Fixed-step','Solver','ode4', ...
    'FixedStep',num2str(Ts),'StopTime',num2str(Tsim));

RMSE = zeros(4,2);  TAUPK = zeros(4,2);  EFIN = zeros(4,2);
E = cell(1,4);  TAU = cell(1,4);  TT = cell(1,4);
for c = 1:4
    in = Simulink.SimulationInput(model);
    in = in.setVariable('q_des_sim',   [t q_des]);
    in = in.setVariable('dq_des_sim',  [t qd_des]);
    in = in.setVariable('ddq_des_sim', [t qdd_des]);
    in = in.setVariable('q0', q0);   in = in.setVariable('dq0', dq0);
    in = in.setVariable('Kp', cfg(c).Kp);   in = in.setVariable('Kd', cfg(c).Kd);
    in = in.setVariable('k_ct', cfg(c).k_ct);
    in = in.setVariable('k_gp', cfg(c).k_gp);
    out = sim(in);
    if ~isempty(out.ErrorMessage), error('%s: %s', names{c}, out.ErrorMessage); end

    tt  = out.tout;
    q   = getsig(out,'q_data');   tau = getsig(out,'tau_data');
    e   = interp1(t,q_des,tt) - q;
    E{c} = e;  TAU{c} = tau;  TT{c} = tt;
    RMSE(c,:)  = sqrt(mean(e.^2));
    TAUPK(c,:) = max(abs(tau));
    EFIN(c,:)  = e(end,:);

    if c == 4      % verifica del GP in anello chiuso
        dq  = getsig(out,'dq_data');  ddq = getsig(out,'ddq_data');
        fgp = getsig(out,'fgp_data');
        tau_hat = zeros(size(q));
        for n = 1:size(q,1)
            tau_hat(n,:) = computed_torque_model(q(n,:).', dq(n,:).', ddq(n,:).', dq(n,:).').';
        end
        tau_tilde = tau - tau_hat;          % tau - tau_hat con stati misurati
        figure; 
        for i = 1:2
            subplot(2,1,i); plot(tt,tau_tilde(:,i),'k',tt,fgp(:,i),'r'); grid on
            ylabel(sprintf('\\tau~_%d (N m)',i)); legend('vero','GP (anello chiuso)');
        end
        xlabel('t (s)');
    end
end

for c = 1:4
    fprintf('\n%s\n', names{c});
    for j = 1:2
        [pk, idx] = max(abs(TAU{c}(:,j)));
        fprintf('Giunto %d: picco %.2f Nm a t = %.3f s\n', ...
            j, pk, TT{c}(idx));
    end
end

%% 4) Risultati
T = table(names', RMSE(:,1), RMSE(:,2), EFIN(:,1), EFIN(:,2), TAUPK(:,1), TAUPK(:,2), ...
    'VariableNames',{'ctrl','RMSE_q1','RMSE_q2','e_fin1','e_fin2','tau_pk1','tau_pk2'});
disp(T);
fprintf('CT-GP vs CT: RMSE ridotto del %.0f %% (q1) e %.0f %% (q2)\n', ...
    100*(1-RMSE(4,1)/RMSE(3,1)), 100*(1-RMSE(4,2)/RMSE(3,2)));


figure; % grafico a barre
b = bar(RMSE);
set(gca, 'XTick', 1:numel(names), ...
         'XTickLabel', names, ...
         'YScale', 'log');
grid on;
ylabel('RMSE (rad)');
legend('Joint 1', 'Joint 2', 'Location', 'best');

% Aggiunge i valori sopra ogni barra
for j = 1:size(RMSE,2)
    for i = 1:size(RMSE,1)
        valore = RMSE(i,j);

        text(b(j).XEndPoints(i), valore*1.15, ...
            sprintf('%.2g', valore), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9);
    end
end

figure;                                       % errori nel tempo
for i = 1:2
    subplot(2,1,i); hold on
    for c = 1:4, plot(TT{c}, E{c}(:,i)); end
    grid on; ylabel(sprintf('e_%d (rad)',i)); legend(names);
end
xlabel('t (s)');

figure; hold on                               % coppia sul giunto 1, come Fig. 8
for c = 1:4, plot(TT{c}, TAU{c}(:,1)); end
grid on; legend(names); xlabel('t (s)'); ylabel('\tau_1 (N m)');

figure; hold on                               % coppia sul giunto 2
for c = 1:4, plot(TT{c}, TAU{c}(:,2)); end
grid on; legend(names); xlabel('t (s)'); ylabel('\tau_2 (N m)');

save('risultati_confronto.mat','T','RMSE','TAUPK','EFIN');

function X = getsig(out,name)
    s = out.get(name);
    if isa(s,'timeseries'), X = s.Data; else, X = s; end
    if ndims(X) == 3, X = squeeze(X).'; end
end