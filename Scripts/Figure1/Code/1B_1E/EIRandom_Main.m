%% EI Random NETWORK (Main)
addpath('./EI_Functions');
addpath('./Helper_Functions');
for Option = 1
    for Iteration = 1:2
        SimInputs = 1; %Simulate Inputs (0 If Loading Old Inputs)
        rng(Iteration)
        
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
        EINetwork;

        % Compute And Plot Results
        EIResults;

        if Iteration==1
            Nidx1=neuron_indices;
            FiringRate1=FiringRate;
        else
            Nidx2=neuron_indices;
            FiringRate2=FiringRate;
        end

    end

%% Remapping Results
figure; 
subplot(1,3,1);imagesc(FiringRate1(Nidx1,n_pos+1:n_pos*6),[0 1]);
title('Context A');xlabel('Time');ylabel('Neuron (Sorting A)')
subplot(1,3,2);imagesc(FiringRate2(Nidx1,n_pos+1:n_pos*6),[0 1])
title('Context B');xlabel('Time');ylabel('Neuron (Sorting A)')
subplot(1,3,3);imagesc(FiringRate2(Nidx2,n_pos+1:n_pos*6),[0 1])
title('Context B');xlabel('Time');ylabel('Neuron (Sorting B)')

end
