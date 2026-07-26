%% EI Random NETWORK (Main)z
clear;close all; 
addpath('./EI_Functions');
addpath('./Helper_Functions');

SlowTm=0;
Save=0;

rng(0)
for Option = 1
    for Iteration = 1
        SimInputs = 1; %Simulate Inputs (0 If Loading Old Inputs)
        
        % Parameters
        n_laps = 5; % Each Lap = n_pos sec
        n_pos = 30;
        n_input = 1000; % Tuned input neurons
        n_excit = 500; % Excitatory neurons
        n_inhib = 50; % Inhibitory neurons
        mean_dt = 0.0005; %Integration Time
        sigma_dt = 0; 
        mean_v = 1; % Running Velocity
        sigma_v = 0;

        %% Initialize LIF
        EIInitialize; %In EI_Functions Folder
        EIPositions; %In EI_Functions Folder

        % Initialize storage matrices
        spike_mat_excit = zeros(n_excit, length(total_time_vec));
        spike_mat_inhib = zeros(n_inhib, length(total_time_vec));
        InputSpikes = zeros(n_input, length(total_time_vec));

        excit_spikes = zeros(1, n_excit);                       % stores the spiking at each time step of the excitatory neurons 
        inhib_spikes = zeros(1, n_inhib);                       
        excit_spikes_2 = zeros(1, n_excit);                     % buffer that stores the t+1 activity so we can still keep the t activity
        inhib_spikes_2 = zeros(1, n_inhib); 
        input_most_recent_fire_times_vec = -100 * ones(1, n_input);
        excit_most_recent_fire_times_vec = -100 * ones(1, n_excit);
        inhib_most_recent_fire_times_vec = -100 * ones(1, n_inhib);

        %% Run The Network
        EINetwork;

        % Compute And Plot Results
        EIResults;

    end
end


if Save==1
%% ============================================================
% Save simulation output for post hoc analysis
%% ============================================================

save_dir = ['/Users/josehurtado/Documents/MATLAB/Final_Manuscript/' ...
    'Random_EIPlaceNet_Internal/Supplementary/S04/Conductance_Model/' ...
    'PreRun_Data/Fast_Pos'];

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save_file = fullfile(save_dir, 'ConductanceModel_Data.mat');

save(save_file, 'total_time_vec', ...
    'dt_vec', ...
    'VmE', ...
    'VmE_pre', ...
    'VinfE', ...
    'tau_eff_E', ...
    'gInputE_history', ...
    'gEE_history', ...
    'gExcE', ...
    'gInhE', ...
    'gTotE', ...
    'IInputE', ...
    'IEE', ...
    'IExcE', ...
    'IInhE', ...
    'ILeakE', ...
    'INetE', ...
    'InputSpikes', ...
    'spike_mat_excit', ...
    'spike_mat_inhib', ...
    'InputRate_t', ...
    'positions', ...
    'pf_cell', ...
    'Ordering', ...
    'W_inputE', ...
    'C', ...
    '-v7.3');
    
fprintf('\nSaved conductance-model data to:\n%s\n', save_file);
end