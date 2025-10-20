%% EI Random NETWORK (Main)
addpath('./EI_Functions');
addpath('./Helper_Functions');

PreRunData=1; 

%% PreRun Colors
%Colors=['g','m','c','y','r']; ID=[50 I, 200 I, 500 I, 1000 I, 1000 I (All same weights)]


if PreRunData==0
%% Set Loops
Data = cell(1, 4); % 1 row, 3 columns (adjust dimensions as needed)
Sparsitydata=cell(1, 4);

for Instance = 1:5

Instance

% AcrossInhibition == 1
% AcrossInputNoise == 2
% AcrossOutputNoise == 3

Across = 1;

%% Begin Loop
Sparsity=zeros(21,1);
data=zeros(size(Sparsity));
Num_I=[50,250,500,1000,1000];
for Iteration = 1:21 % Main Loop

    Iteration

for Option = 1:10
        SimInputs = 1; %Simulate Inputs (0 If Loading Old Inputs)
        
        % Parameters
        n_laps = 5; % Each Lap = n_pos sec
        n_pos = 30;
        n_input = 1000; % Tuned input neurons
        n_excit = 1000; % Excitatory neurons
        n_inhib = Num_I(Instance); % Inhibitory neurons
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

       % Compute the proportion of rows that are nonzero for each column (time point)
        nonzero_proportion = sum(FiringRate ~= 0, 1) / size(FiringRate, 1);
        
       % Compute the average of these proportions over all time points
        Sparsity(Iteration,Option) = 1-mean(nonzero_proportion);
        
        Run_Decoder1Fi; %In Helper_Functions Folder




end
end

Data{Instance} = data;
Sparsitydata{Instance} = Sparsity;


end

elseif PreRunData==1
    addpath('./Data');
    addpath('./EI_Functions');
    addpath('./Helper_Functions');
    
    load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S01/Data/Data.mat')
    load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S01/Data/Sparsity.mat')

end

for Instance = 1:5


data = Data{Instance};
EI_FullAnalysis_Results; %In EI_Functions Folder
end

for Instance = 1:5
data = Sparsitydata{Instance};
EI_FullAnalysis_Results; %In EI_Functions Folder
end

