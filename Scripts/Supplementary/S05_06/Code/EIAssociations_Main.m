addpath('./EI_Functions');
addpath('./Helper_Functions');

%% EIAssociations_Main - Runs Training/Tagging, And Recall Simulations

numSims=10; %Number of simulations

Reversals = zeros(numSims,3);
Freezes = zeros(numSims,3);
SpecialPositions = zeros(numSims,1);
TaggedSize = zeros(numSims,1);
CollectIndices= cell(2,numSims);
CollectRates= cell(1,numSims);
CollectReverse = cell(1, numSims);
CollectFreeze = cell(1, numSims);
CollectAvoidRates = cell(1,numSims);

Recall=1; 

% Specify the indices of the files you want to use (choose any 10 out of 12)
selectedIndices = 1:numSims;  % Choose 10 Separate Simulations
FigureToPlot = 1;

Sims=1:numSims;

for J = Sims

    J
    
    %% Parameters 
    StartGoal = 0; %Goal will be in second track
    
    EIInitialize;
    EIPositions;
    
    
    %% Initializing matrices & vectors to store various information about the dynamic
    spike_mat_excitf = zeros(n_excit, length(total_time_vec)); % stores excitatory spikes 
    spike_mat_inhibf = zeros(n_inhib, length(total_time_vec));
    inhib_cumul_inputf = zeros(n_inhib, length(total_time_vec));
    excit_cumul_inputf = zeros(n_excit, length(total_time_vec));
    NetCurrentET = zeros(n_excit, length(total_time_vec));
    NetCurrentIT = zeros(n_inhib, length(total_time_vec));
    FrustrationET = zeros(n_excit, length(total_time_vec));
    FrustrationIT = zeros(n_inhib, length(total_time_vec));
    excit_spikesf = zeros(1, n_excit);                       % stores the spiking at each time step of the excitatory neurons 
    inhib_spikesf = zeros(1, n_inhib);                       
    excit_spikes_2f = zeros(1, n_excit);                     % buffer that stores the t+1 activity so we can still keep the t activity
    inhib_spikes_2f = zeros(1, n_inhib); 
    excit_cum_inputf = zeros(1, n_excit);                    % stores the voltage of the excitatory neurons
    inhib_cum_inputf = zeros(1, n_inhib); 
    x_excitf= zeros(1, n_excit);                            % sotres the exponentials used in calculating the STDP changes
    x_inhibf = zeros(1, n_inhib);
    input_most_recent_fire_times_vec = zeros(1, n_input) - 100; % stores the most recent firing times of the input neurons 
    excit_most_recent_fire_times_vecf = zeros(1, n_excit) - 100;
    inhib_most_recent_fire_times_vecf = zeros(1, n_inhib) - 100;
    AveFrustIaT = zeros(1,length(total_time_vec));
    AveFrustIbT = zeros(1,length(total_time_vec));
    AveFrustET = zeros(1,length(total_time_vec));
    Input_Spikes = zeros(length(total_time_vec),n_input);
    
    pf_cell=pf_cell1;
    
    EINetwork;
    EIResultsandRecall;



% Collect data
CollectReverse{J} = ReverseData;
CollectFreeze{J} = FreezeData;
CollectPositions(:, J) = PosRecall;

end



SummaryPlots;

