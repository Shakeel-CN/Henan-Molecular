% Define parameters
k_prod = 1;                % ABA production rate constant (concentration/time)
k_deg = 0.1;               % ABA degradation rate constant (1/time)
k_on = 1e4;                % ABA binding to receptor, association rate constant (M^-1 s^-1)
k_off = 1e-3;              % ABA dissociating from receptor, dissociation rate constant (s^-1)
k_half = 1e-2;             % Half-maximal concentration for the Hill equation
n = 0.4;                   % Cooperativity coefficient in ligand binding for Hill equation
k_trans = 1;               % Transcription rate constant for sensor protein mRNA (concentration/time)
k_trans_deg = 0.1;         % Degradation rate constant for sensor protein mRNA (1/time)
k_trans_syn = 1;           % Translation rate constant for sensor protein synthesis (concentration/time)
k_deg_prot = 0.1;          % Degradation rate constant for sensor protein (1/time)
k_response = 0.5;          % Rate constant for bioluminescence response (1/time)
R_total = 1;               % Total receptor concentration (concentration)

% Define initial conditions
C_ABA_init = 0:5:50;       % Initial ABA concentration range (0 to 50 µM)
C_ABA_R_init = 0;          % Initial concentration of ABA-receptor complex
mRNA_init = 0;             % Initial concentration of sensor protein mRNA
Protein_init = 0;          % Initial concentration of sensor protein
I_Bioluminescence_init = 0; % Initial bioluminescence response

% Time span for the simulation
tspan = [0 100];

% Define parameter struct
params = struct('k_prod', k_prod, 'k_deg', k_deg, 'k_on', k_on, 'k_off', k_off, ...
                'k_half', k_half, 'n', n, 'k_trans', k_trans, 'k_trans_deg', k_trans_deg, ...
                'k_trans_syn', k_trans_syn, 'k_deg_prot', k_deg_prot, 'k_response', k_response, ...
                'R_total', R_total);

% Preallocate results for final bioluminescence response and detection time
final_bioluminescence = zeros(length(C_ABA_init), 1);
detection_time = zeros(length(C_ABA_init), 1);

% Define a threshold for bioluminescence detection
threshold = 0.5;

% Solve the ODEs for different initial ABA concentrations
for i = 1:length(C_ABA_init)
    % Initial conditions for the current simulation
    Y0 = [C_ABA_init(i), C_ABA_R_init, mRNA_init, Protein_init, I_Bioluminescence_init];
    
    % Solve ODEs
    [T, Y] = ode15s(@(t, Y) aba_biosensor_odes(t, Y, params), tspan, Y0);
    
    % Store final bioluminescence response
    final_bioluminescence(i) = Y(end, 5);
    
    % Find the detection time when bioluminescence exceeds the threshold
    idx = find(Y(:, 5) >= threshold, 1);
    if ~isempty(idx)
        detection_time(i) = T(idx);
    else
        detection_time(i) = NaN; % If the threshold is not reached
    end
end

% Plot final bioluminescence response
figure;
plot(C_ABA_init, final_bioluminescence, '-o');
xlabel('Initial ABA Concentration (µM)');
ylabel('Final Bioluminescence Response');
title('Bioluminescence Response for Different Initial ABA Concentrations');
grid on;

% Plot detection time
figure;
plot(C_ABA_init, detection_time, '-o');
xlabel('Initial ABA Concentration (µM)');
ylabel('Detection Time (s)');
title('Detection Time for Different Initial ABA Concentrations');
grid on;

% Define the system of ODEs
function dYdt = aba_biosensor_odes(t, Y, params)
    % Unpack parameters
    k_prod = params.k_prod;
    k_deg = params.k_deg;
    k_on = params.k_on;
    k_off = params.k_off;
    k_half = params.k_half;
    n = params.n;
    k_trans = params.k_trans;
    k_trans_deg = params.k_trans_deg;
    k_trans_syn = params.k_trans_syn;
    k_deg_prot = params.k_deg_prot;
    k_response = params.k_response;
    R_total = params.R_total;
    
    % Unpack state variables
    C_ABA = Y(1);
    C_ABA_R = Y(2);
    mRNA = Y(3);
    Protein = Y(4);
    I_Bioluminescence = Y(5);
    
    % Define ODEs
    dC_ABA_dt = k_prod - k_deg * C_ABA - k_on * C_ABA * Protein + k_off * C_ABA_R;
    dC_ABA_R_dt = k_on * C_ABA * Protein - k_off * C_ABA_R;
    dmRNA_dt = k_trans * (C_ABA^n / (k_half^n + C_ABA^n)) - k_trans_deg * mRNA;
    dProtein_dt = k_trans_syn * mRNA - k_deg_prot * Protein;
    dI_Bioluminescence_dt = k_response * (C_ABA_R^n / (k_half^n + C_ABA_R^n));
    
    % Pack derivatives
    dYdt = [dC_ABA_dt; dC_ABA_R_dt; dmRNA_dt; dProtein_dt; dI_Bioluminescence_dt];
end

