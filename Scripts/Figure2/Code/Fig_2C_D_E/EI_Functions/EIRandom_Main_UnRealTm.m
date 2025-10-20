%% EI Random NETWORK (Main)
addpath('./EI_Functions');
addpath('./Helper_Functions');

%% Set Loops
SlowTm = 1;
N_params = 21;
N_samples = 1;
X1=cell(N_params,N_samples);

for Instance = 1


% AcrossInhibition == 1
% AcrossInputNoise == 2
% AcrossOutputNoise == 3

Across = Instance; % Select Loop (Across inhibition or noise)

%% Begin Loop
Sparsity=zeros(N_params,1);
data=zeros(size(Sparsity));

for Iteration = 1:N_params % Main Loop

    Iteration

for Option = 1:N_samples
        SimInputs = 1; %Simulate Inputs (0 If Loading Old Inputs)
        
        % Parameters
        n_laps = 5; % Each Lap = n_pos sec
        n_pos = 30;
        n_input = 1000; % Tuned input neurons
        n_excit = 500; % Excitatory neurons
        n_inhib = 50; % Inhibitory neurons
        mean_dt = 0.005; %Integration Time
        sigma_dt = 0; 
        mean_v = 1; % Running Velocity
        sigma_v = 0;

        %% Initialize LIF
        
        EIInitialize; %In EI_Functions Folder
        EIPositions; %In EI_Functions Folder

        % Initialize storage matrices
        spike_mat_excit = zeros(n_excit, length(total_time_vec));
        spike_mat_inhib = zeros(n_inhib, length(total_time_vec));
        VmE = zeros(n_excit, length(total_time_vec));
        VmI = zeros(n_inhib, length(total_time_vec));
        InputSpikes = zeros(n_input, length(total_time_vec));
        excit_cum_input = zeros(1, n_excit);
        inhib_cum_input = zeros(1, n_inhib);
        excit_spikes = zeros(1, n_excit);                       % stores the spiking at each time step of the excitatory neurons 
        inhib_spikes = zeros(1, n_inhib);                       
        excit_spikes_2 = zeros(1, n_excit);                     % buffer that stores the t+1 activity so we can still keep the t activity
        inhib_spikes_2 = zeros(1, n_inhib); 
        input_most_recent_fire_times_vec = -100 * ones(1, n_input);
        excit_most_recent_fire_times_vec = -100 * ones(1, n_excit);
        inhib_most_recent_fire_times_vec = -100 * ones(1, n_inhib);
        

        %% Run The Network
        EINetwork; %In EI_Functions Folder

        % Compute And Plot Results
        EIResults; %In EI_Functions Folder


X1{Iteration,Option}=FiringRate;

end


end

end

