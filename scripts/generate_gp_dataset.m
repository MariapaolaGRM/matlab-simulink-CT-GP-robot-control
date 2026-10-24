%% Genera il dataset di addestramento e validazione per il Gaussian Process

%clear; close all; clc; 

root = fileparts(fileparts(mfilename('fullpath')));   % cartella che contiene Scripts, Simulink, src
addpath(fullfile(root, 'src'), fullfile(root, 'Simulink'), fullfile(root, 'data'));
%% Impostazioni
model   = 'GP_DataCollection_2R';     % nome modello di training
Ts      = 1e-3;                   % passo del solver (1ms) [s]
Ts_data = 30e-3;                  % passo di registrazione dei dati [s]
t_tail  = 0.5;                    % tempo extra dopo la fine del moto [s]
q_lim   = [0 pi; -pi/4 3*pi/4];   % [min max] di q1 e q2: deve contenere la zona di lavoro
T_range = [1.5 5];                % durate dei movimenti [s]
N_train = 40;                     % num traiettorie di training
frac_quintic = 0.6;               % 60% dei movimenti sono quintici (40% trapezoidali)
N_val   = 6;                      % num traiettorie di validazione

%% Guadagni HG-PD 
p_hat = get_params(false);        % parametri del modello stimato
q2s = linspace(-pi, pi, 200);     % genera 200 valori di q2
Hd  = zeros(2, numel(q2s));       % per ogni valore valuta la diagonale della matr di inerzia stimata
for k = 1:numel(q2s)
    Hd(:,k) = diag(H([0; q2s(k)], p_hat));
end
H_nom = min(Hd, [], 2);           % seleziona il min di ogni elemento diagonale
wn = 20;  
delta = 1;
Kp = diag(H_nom * wn^2);
Kd = diag(2*delta*wn*H_nom);

%% Solver
load_system(model);
set_param(model,'SolverType','Fixed-step','Solver','ode4','FixedStep',num2str(Ts));

%% Training
rng(1); % fissa il generatore casuale per rendere il dataset riproducibile
nq = round(frac_quintic*N_train);  % num traiett quintiche
types_train = [repmat({'quintic'},1,nq), repmat({'trapezoid'},1,N_train-nq)];
D_train = collect_dataset(model, types_train, q_lim, T_range, Ts, Ts_data, t_tail, Kp, Kd);

%% Validazione
rng(2);
D_val = collect_dataset(model, repmat({'cubic'},1,N_val), q_lim, T_range, Ts, Ts_data, t_tail, Kp, Kd);

%% Salvataggio dataset
dataDir = fullfile(root, 'data');
%if ~exist(dataDir, 'dir'), mkdir(dataDir); end          % crea la cartella se non esiste

save(fullfile(dataDir, 'dataset.mat'), 'D_train', 'D_val','Kp','Kd','Ts_data','q_lim');
fprintf('Training: %d punti | Validazione: %d punti\n', size(D_train.X,1), size(D_val.X,1));

%% Controllo di copertura dello spazio degli ingressi
X = D_train.X;      % colonne: [ddq1 ddq2 dq1 dq2 q1 q2]
figure;
subplot(1,3,1); plot(X(:,5),X(:,6),'.'); grid on; xlabel('q_1'); ylabel('q_2'); title('configurazioni')
subplot(1,3,2); plot(X(:,3),X(:,1),'.'); grid on; xlabel('dq_1'); ylabel('ddq_1'); title('giunto 1')
subplot(1,3,3); plot(X(:,4),X(:,2),'.'); grid on; xlabel('dq_2'); ylabel('ddq_2'); title('giunto 2')

%% FUNZIONI 
function D = collect_dataset(model, types, q_lim, T_range, Ts, Ts_data, t_tail, Kp, Kd)
    nT   = numel(types);         % num di elementi in types
    step = round(Ts_data/Ts);
    X = [];  Y = [];  TAU = [];  ID = [];  TIME = [];
    for k = 1:nT
        T   = T_range(1) + diff(T_range)*rand;       % durata della traiettoria estratta dal range
        tcf = 0.25 + 0.15*rand;                      % frazione di accelerazione (trapezoidale) tra 25% e 40%
        [qi, qf] = random_waypoints(q_lim);
        Tsim = T + t_tail;                           % tempo di simulazione (comprende fase cost finale)
        t = (0:Ts:Tsim)';
        [q_d, dq_d, ddq_d] = traj_eval(types{k}, qi, qf, T, t, tcf);

        in = Simulink.SimulationInput(model);
        in = in.setVariable('q_des_sim',   [t q_d]);
        in = in.setVariable('dq_des_sim',  [t dq_d]);
        in = in.setVariable('ddq_des_sim', [t ddq_d]);
        in = in.setVariable('q0',  qi(:));           % condizione iniziale integratore di q
        in = in.setVariable('dq0', [0; 0]);          % condizione iniziale integratore di dq
        in = in.setVariable('Kp', Kp);
        in = in.setVariable('Kd', Kd);
        in = in.setModelParameter('StopTime', num2str(Tsim));

        out = sim(in);
        if ~isempty(out.ErrorMessage)
            warning('Traiettoria %d fallita: %s', k, out.ErrorMessage);  continue
        end

        q   = getsig(out,'q_data');
        dq  = getsig(out,'dq_data');
        ddq = getsig(out,'ddq_data');
        tt  = getsig(out,'tau_tilde_data');          % tau - tau_hat
        tau = getsig(out,'tau_data');                % coppia applicata
        idx = 1:step:size(q,1);

        X   = [X;   ddq(idx,:) dq(idx,:) q(idx,:)];  
        Y   = [Y;   tt(idx,:)];                      
        TAU = [TAU; tau(idx,:)];                     
        ID  = [ID;  k*ones(numel(idx),1)];           
        TIME= [TIME; out.tout(idx)];                 
        fprintf('%2d/%d  %-9s T=%.2f s  q:[%.2f %.2f]->[%.2f %.2f]\n', ...
                k, nT, types{k}, T, qi, qf);
    end
    D = struct('X',X,'Y',Y,'tau',TAU,'id',ID,'t',TIME,'types',{types});
end

function [qi, qf] = random_waypoints(q_lim)
    lo = q_lim(:,1);  % limite inferiore 
    hi = q_lim(:,2);  % limite superiore
    qi = lo + (hi-lo).*rand(2,1); % estrae config iniziale casualmente
    qf = lo + (hi-lo).*rand(2,1); % anche qf casuale
    
    qf = min(max(qf,lo),hi);
end

function [q, v, a] = traj_eval(type, qi, qf, T, t, tcf)
    ti = 0;  tf = T;
    
    n = numel(t);  
    q = zeros(n,2);  
    v = q;  
    a = q;
    tt = min(t, T);         % dopo T il riferimento resta fermo
    c  = ones(size(tt));
    for j = 1:2
        qi_j = qi(j);  qf_j = qf(j);  D = qf_j - qi_j;
        vi = 0;  vf = 0;  aci = 0;  acf = 0;
        switch type
            case 'quintic'
                Tm = [1 ti ti^2 ti^3 ti^4 ti^5; 0 1 2*ti 3*ti^2 4*ti^3 5*ti^4; ...
                      0 0 2 6*ti 12*ti^2 20*ti^3; ...
                      1 tf tf^2 tf^3 tf^4 tf^5; 0 1 2*tf 3*tf^2 4*tf^3 5*tf^4; ...
                      0 0 2 6*tf 12*tf^2 20*tf^3];
                b  = [qi_j; vi; aci; qf_j; vf; acf];
                p  = Tm\b;
                q(:,j) = p(1)*c + p(2)*tt + p(3)*tt.^2 + p(4)*tt.^3 + p(5)*tt.^4 + p(6)*tt.^5;
                v(:,j) = p(2)*c + 2*p(3)*tt + 3*p(4)*tt.^2 + 4*p(5)*tt.^3 + 5*p(6)*tt.^4;
                a(:,j) = 2*p(3)*c + 6*p(4)*tt + 12*p(5)*tt.^2 + 20*p(6)*tt.^3;
            case 'cubic'
                Tm = [1 ti ti^2 ti^3; 0 1 2*ti 3*ti^2; ...
                      1 tf tf^2 tf^3; 0 1 2*tf 3*tf^2];
                b  = [qi_j; vi; qf_j; vf];
                p  = Tm\b;
                q(:,j) = p(1)*c + p(2)*tt + p(3)*tt.^2 + p(4)*tt.^3;
                v(:,j) = p(2)*c + 2*p(3)*tt + 3*p(4)*tt.^2;
                a(:,j) = 2*p(3)*c + 6*p(4)*tt;
            case 'trapezoid'
                tc = tcf*T;  vm = D/(T - tc);  ac = vm/tc;   % con segno
                i1 = tt < tc;  i3 = tt > T - tc;  i2 = ~i1 & ~i3;
                q(i1,j) = qi(j) + 0.5*ac*tt(i1).^2;
                q(i2,j) = qi(j) + vm*(tt(i2) - tc/2);
                q(i3,j) = qf(j) - 0.5*ac*(T - tt(i3)).^2;
                v(i1,j) = ac*tt(i1);
                v(i2,j) = vm;
                v(i3,j) = ac*(T - tt(i3));
                a(i1,j) =  ac;
                a(i3,j) = -ac;
            otherwise
                error('Tipo di traiettoria sconosciuto: %s', type);
        end
        a(:,j) = a(:,j) .* (t < T);   % accelerazione nulla a fine moto
    end
end

function X = getsig(out, name)
    s = out.get(name);
    if isa(s,'timeseries')
        X = s.Data; 
    else 
        X = s; 
    end
    if ndims(X) == 3
        X = squeeze(X).'; 
    end 
end