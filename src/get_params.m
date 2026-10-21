function p = get_params(is_real)
    % --- Parametri nominali ---
    a1 = 1;   a2 = 1;       % lunghezza link [m]
    g0 = 9.81;              % accelerazione gravitazionale [m/s^2]
    l1 = a1/2;  l2 = a2/2;  % posizione baricentro link [m]
    ml1 = 50;  ml2 = 50;    % massa link [kg]
    mm1 = 5;   mm2 = 5;     % massa motore [kg]
    Il1 = 10;  Il2 = 10;    % inerzia link [kg m^2]
    Im1 = 0.01; Im2 = 0.01; % inerzia rotori [kg m^2]
    kr1 = 100; kr2 = 100;   % rapporti di riduzione
    Fv = [0; 0];            % coefficienti di attrito viscoso
    % Fc = [0; 0];

    if is_real
        % comprende dinamiche non modellate del robot reale
        Fv = [2; 2];
        % Fc = [1; 1];
        % mp  = 0; % nessun payload 5;           % payload
        % m2t = ml2 + mp;
        % l2e = (ml2*l2 + mp*a2) / m2t;
        % Il2 = Il2 + ml2*(l2e - l2)^2 + mp*(a2 - l2e)^2;
        % ml2 = m2t;
        % l2  = l2e;
    else
        % modello stimato: errore parametrico, niente attrito 
        e = 0.15;
        ml1 = ml1*(1+e);   ml2 = ml2*(1-e);
        mm1 = mm1*(1+e);   mm2 = mm2*(1-e);
        Il1 = Il1*(1+e);   Il2 = Il2*(1-e);
    end

    % --- Costruzione della struct in un colpo solo ---
    p = struct( ...
        'a1',  a1,  'a2',  a2,  'g',   g0, ...
        'l1',  l1,  'l2',  l2, ...
        'ml1', ml1, 'ml2', ml2, ...
        'mm1', mm1, 'mm2', mm2, ...
        'Il1', Il1, 'Il2', Il2, ...
        'Im1', Im1, 'Im2', Im2, ...
        'kr1', kr1, 'kr2', kr2, ...
        'Fv',  Fv ); %,  'Fc',  Fc ... );
end