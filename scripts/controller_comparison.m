% Confronta le prestazioni di diverse strategie di controllo del manipolatore
% robotico a due giunti mediante simulazioni in MATLAB/Simulink.

%close all; clc;

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'));

model = 'CT_PD_robot_2R';
Ts = 1e-3;  tf = 3;  Tsim = 6;
%% Guadagni LG e HG
p_hat = get_params(false);
q2s = linspace(-pi,pi,200);  Hd = zeros(2,numel(q2s));
for k = 1:numel(q2s), Hd(:,k) = diag(H([0;q2s(k)],p_hat)); end
H_nom = min(Hd,[],2);  zeta = 1;
Kp_LG = diag(H_nom*4^2);   Kd_LG = diag(2*zeta*4*H_nom);
Kp_HG = diag(H_nom*20^2);  Kd_HG = diag(2*zeta*20*H_nom);

%% Traiettoria di test (cubica, tenuta ferma dopo tf)
q0 = [pi/4; pi/8];  qf = [pi/2; pi/4];  dq0 = [0;0];
t  = (0:Ts:Tsim)';  tc = min(t,tf);  on = double(t <= tf);
q_des = zeros(numel(t),2);  qd_des = q_des;  qdd_des = q_des;
for j = 1:2
    D = qf(j)-q0(j);  a2 = 3*D/tf^2;  a3 = -2*D/tf^3;
    q_des(:,j)   = q0(j) + a2*tc.^2 + a3*tc.^3;
    qd_des(:,j)  = (2*a2*t + 3*a3*t.^2).*on;
    qdd_des(:,j) = (2*a2 + 6*a3*t).*on;
end

% Visualizzazione della traiettoria desiderata
% figure;
% for j = 1:2
%     subplot(3,1,1);
%     plot(t,q_des(:,j)); hold on; grid on;
%     ylabel('q_d (rad)');
% 
%     subplot(3,1,2);
%     plot(t,qd_des(:,j)); hold on; grid on;
%     ylabel('dq_d (rad/s)');
% 
%     subplot(3,1,3);
%     plot(t,qdd_des(:,j)); hold on; grid on;
%     ylabel('ddq_d (rad/s^2)');
%     xlabel('t (s)');
% end
%% Quattro controllori
% Configurazione del modello Simulink
load_system(model);
set_param(model, ...
    'SolverType', 'Fixed-step', ...
    'Solver', 'ode4', ...
    'FixedStep', num2str(Ts), ...
    'StopTime', num2str(Tsim));

% Configurazione dei controllori
names = {'LG-PD','HG-PD','CT','CT-GP'};

cfg = [ ...
    struct('Kp', Kp_LG, 'Kd', Kd_LG, 'k_ct', 0, 'k_gp', 0), ...
    struct('Kp', Kp_HG, 'Kd', Kd_HG, 'k_ct', 0, 'k_gp', 0), ...
    struct('Kp', Kp_LG, 'Kd', Kd_LG, 'k_ct', 1, 'k_gp', 0), ...
    struct('Kp', Kp_LG, 'Kd', Kd_LG, 'k_ct', 1, 'k_gp', 1)];

% Definizione degli scenari
scen = struct('nome',{'ideale','disturbo','rumore'}, 'kd',{0,1,0}, 'kn',{0,0,1});

nseed = 10;    % ripetizioni solo per il rumore
nscen = numel(scen);
nctrl = numel(names);

% Massimo numero di ripetizioni previsto
nrep_max = nseed;

% Metriche: scenario x controllore x ripetizione x giunto
RMSE  = nan(nscen, nctrl, nrep_max, 2);
EFIN  = nan(nscen, nctrl, nrep_max, 2);
TAUPK = nan(nscen, nctrl, nrep_max, 2);
TAURMS = nan(nscen, nctrl, nrep_max, 2);

% Risultati completi delle simulazioni
E   = cell(nscen, nctrl, nrep_max);
TAU = cell(nscen, nctrl, nrep_max);
TT  = cell(nscen, nctrl, nrep_max);

% Simulazioni
for s = 1:nscen
  for c = 1:nctrl                      % LG-PD, HG-PD, CT, CT-GP
    if scen(s).kn == 1
        nrep = nseed;
    else
        nrep = 1;
    end
    for r = 1:nrep
      in = Simulink.SimulationInput(model);
      in = in.setVariable('q_des_sim',   [t q_des]);
      in = in.setVariable('dq_des_sim',  [t qd_des]);
      in = in.setVariable('ddq_des_sim', [t qdd_des]);
     
      in = in.setVariable('q0', q0);% Condizioni iniziali q
      in = in.setVariable('dq0', dq0);% Condizioni iniziali dq

      in = in.setVariable('Kp', cfg(c).Kp);
      in = in.setVariable('Kd', cfg(c).Kd);
      in = in.setVariable('k_ct', cfg(c).k_ct);
      in = in.setVariable('k_gp', cfg(c).k_gp);

      in = in.setVariable('k_dist', scen(s).kd);
      in = in.setVariable('k_noise', scen(s).kn);

      in = in.setVariable('seed', 100*r); % seed del rumore

      out = sim(in);
      if ~isempty(out.ErrorMessage)
          error('Scenario %s, controllore %s: %s', ...
              scen(s).nome, names{c}, out.ErrorMessage);
      end

      tt  = out.tout; % segnali di output
      q   = getsig(out, 'q_data');
      tau = getsig(out, 'tau_data');

      if s == 1 && c == 4 && r == 1    % CT-GP, scenario ideale
          q_gp = q;  tau_gp = tau;  tt_gp = tt;
          dq_gp  = getsig(out,'dq_data');
          ddq_gp = getsig(out,'ddq_data');
          fgp_gp = getsig(out,'fgp_data');
      end
      
      e   = interp1(t,q_des,tt) - q; % errore di inseguimento sulla posizione reale

      % Salvataggio dei segnali
      E{s,c,r}   = e;
      TAU{s,c,r} = tau;
      TT{s,c,r}  = tt;

      % Metriche per ciascun giunto
      RMSE(s,c,r,:) = reshape(sqrt(mean(e.^2, 1)), 1,1,1,2);
      EFIN(s,c,r,:) = reshape(e(end,:), 1,1,1,2);
      TAUPK(s,c,r,:) = reshape(max(abs(tau), [], 1), 1,1,1,2);
      TAURMS(s,c,r,:) = reshape(sqrt(mean(tau.^2, 1)), 1,1,1,2);

      fprintf('Scenario: %s | Controllore: %s | Ripetizione: %d\n', ...
            scen(s).nome, names{c}, r);
    end
  end
end

tau_hat = zeros(size(q_gp));
for n = 1:size(q_gp,1)
    tau_hat(n,:) = computed_torque_model(q_gp(n,:).', dq_gp(n,:).', ddq_gp(n,:).', dq_gp(n,:).').';
end
tau_tilde = tau_gp - tau_hat;                  % residuo vero, con gli stati misurati
figure('Name','Residuo vero e GP in anello chiuso');
for i = 1:2
    subplot(2,1,i); plot(tt_gp, tau_tilde(:,i),'k', tt_gp, fgp_gp(:,i),'r'); grid on
    ylabel(sprintf('\\tau~_%d (N m)',i)); legend('vero','GP (anello chiuso)');
end
xlabel('t (s)');
%% Risultati
for s = 1:nscen
    if scen(s).kn == 1 % num ripetizioni per scenario
        nrep = nseed;
    else
        nrep = 1;
    end
    Rm = reshape(mean(RMSE(s,:,1:nrep,:), 3), nctrl, 2);
    Rs = reshape(std(RMSE(s,:,1:nrep,:), 0, 3), nctrl, 2);
    Em = reshape(mean(EFIN(s,:,1:nrep,:), 3), nctrl, 2);
    Pm = reshape(mean(TAUPK(s,:,1:nrep,:), 3), nctrl, 2);
    Qm = reshape(mean(TAURMS(s,:,1:nrep,:), 3), nctrl, 2);
    
    T = table(names', Rm(:,1), Rm(:,2), Em(:,1), Em(:,2), Pm(:,1), Pm(:,2), Qm(:,1), Qm(:,2), ...
        'VariableNames',{'ctrl','RMSE_q1','RMSE_q2','e_fin1','e_fin2','tau_pk1','tau_pk2','tau_rms1','tau_rms2'});
    fprintf('\n=== Scenario: %s (%d ripetizioni) ===\n', scen(s).nome, nrep);
    disp(T);

    figure('Name',['RMSE - ' scen(s).nome]);
    b = bar(Rm);
    set(gca,'XTickLabel',names,'YScale','log'); grid on;
    ylabel('RMSE (rad)'); legend('Joint 1','Joint 2','Location','best');
    for j = 1:2
        for i = 1:nctrl
            text(b(j).XEndPoints(i), Rm(i,j)*1.15, sprintf('%.2g',Rm(i,j)), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
        end
    end

    % Deviazione standard della RMSE
    if nrep > 1
        Tstd = table(names', Rs(:,1), Rs(:,2), ...
            'VariableNames', ...
            {'ctrl','std_RMSE_q1','std_RMSE_q2'});
        disp(Tstd);
    end

    figure('Name',['Errori - ' scen(s).nome]);   % errori nel tempo (prima ripetizione)
    for i = 1:2
        subplot(2,1,i); hold on
        for c = 1:nctrl, plot(TT{s,c,1}, E{s,c,1}(:,i)); end
        grid on; 
        ylabel(sprintf('e_%d (rad)',i)); 
        legend(names);
    end
    xlabel('t (s)');
    
    figure('Name',sprintf('Coppia - %s',scen(s).nome));
    for j = 1:2    % coppie
        subplot(1,2,j); hold on
        for c = 1:nctrl
            plot(TT{s,c,1}, TAU{s,c,1}(:,j)); 
        end
        grid on; 
        legend(names); 
        xlabel('t (s)'); 
        ylabel(sprintf('\\tau_%d (N m)',j));
    end
end

save('risultati_confronto.mat','RMSE','EFIN','TAUPK','TAURMS','scen','names');

%% FUNZIONI
function X = getsig(out,name)
    s = out.get(name);
    if isa(s,'timeseries'), X = s.Data; else, X = s; end
    if ndims(X) == 3, X = squeeze(X).'; end
end