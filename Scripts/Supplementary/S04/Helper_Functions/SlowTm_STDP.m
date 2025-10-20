addpath('./Helper_Functions');
addpath('./SavedFiles');

%% Adjust Save Directory Below (Lines 280-308 and Lines 493-517)

rng(0)

%% Variables to predefine outputs (If changed, make sure to also changes within parfor!)
NumIterations=1; 
NumRepetitions=1; 
NumNeurons=500;
NumInputs=1000;
NumPos = 30;
NumLaps=50;
ProbPF = 0.2;
InitialWeight=0.005;

pf_cell1 = GET_TRACK_PFs(NumPos,ProbPF, NumInputs);

pf_cell2 = GET_TRACK_PFs(NumPos,ProbPF, NumInputs);

n_laps = 50; 
n_pos = 30;

MaxWeight=0.5;

%% Adjusting Input, IE and EI resistance to current
UnitAdjInput = .00015; 
UnitAdjInhibEI = .075; 
UnitAdjInhibIE = .075/20;
 
%% Parameters 
n_input = 1000;                             % number of tuned input neurons                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
n_excit = 500;                              % number of excitatory neurons
n_inhib = 50;                               % number of inhibititory neurons 

%% Parameters 

prob_E_recurrent_connectivity = 0.25;       % probability of any given excitatory neuron will connect to another excitatory neuron
prob_I_E_connectivity = 0.30;               % probability of any inhibitory neuron connecting to an excitatory neuron
prob_E_I_connectivity = 0.175;
prob_I_I_connectivity = .50;              % probability of any I neuron connecting to an inhibitory neuron                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ; %probability of any excitatory neuron connecting to an inhibitory neuron
initial_weight_max_input = 0.05;            % initial weight maximum input to E
initial_weight_max_EE = 0.05 / prob_E_recurrent_connectivity; % initial weight maximum for E to E
initial_weight_max_EI = 0.05 / prob_I_E_connectivity;  % initial weight maximum for E to I
initial_weight_max_IE = 0.05 / prob_E_I_connectivity; % initial weight maximum for I to E
initial_weight_max_II = 0.05 / prob_I_I_connectivity; % initial weight maximum for I to E

I_F_threshE = -55;                           % integrate and fire threshold 
I_F_threshI = -55;                           % integrate and fire threshold 

VL = -70;
VrE = -65;
VrI = -65;
VLeakE = VL*ones(1,n_excit); %Leak Potential
VLeakI = VL*ones(1,n_inhib); %Leak Potential

%RC Circuit Parameters (To Preserve mV and s Units)
C=.0000005; %Capacitance u(mF) where mF = (Amp/mV) | Note S.I. F = (Amp/V)
gL=.000025/50; %Conductance u(mS) where mS = (C/mV)   | Note S.I. S = (C/V)

Tm=(C/gL); %Resulting units are in seconds


n_steps = n_laps * n_pos;                   % number of steps on the track 
input_refract_length = 0;                   % input neuron refractory length 
excit_refract_length = 0.002;                % excit refractory length
inhib_refract_length = 0.002;                % inhibitory refractory length
W_upper_limit = MaxWeight;                        % maximum weight we allow the E to I weights to take
alpha_P = 0;                             % pink noise strength
pf_rate = 12;                                % input neuron rate (in Hz) at center of pf
pf_width = 10;                              % controls the Gaussian pf width (note it's actually the inverse) 
corr_window_size = 0.25;                    % size of window used for smoothing activity for correlation analysis
dt=0.005;
Tsyn=1;

EI_plast=1;
II_plast=0;
IE_plast=1;
EE_plast=0;
eta=0.0005;

%% Initiliazing the input neurons 

%% Recurrent Weights %%
% Note: W_inputE is not plastic
W_inputE = rand(n_input, n_excit) * initial_weight_max_input;

W_EE = 0;%rand(n_excit, n_excit); 
W_II = 0;%rand(n_inhib, n_inhib); 
W_EI = rand(n_excit, n_inhib); 
W_IE = rand(n_inhib, n_excit);

% W_EE(W_EE > prob_E_recurrent_connectivity) = 0;
% W_EE(diag(ones(1,n_excit)) == 1) = 0;       % no neurons can connect recurrently to themselves
% W_EE=W_EE/max(W_EE,[],'all');
% W_EE = W_EE * initial_weight_max_EE;
% 
% W_II(W_II > prob_I_I_connectivity) = 0;
% W_II(diag(ones(1,n_inhib)) == 1) = 0;       % no neurons can connect recurrently to themselves
% W_II=W_II/max(W_II,[],'all');
% W_II = W_II * initial_weight_max_II;


% W_EI(W_EI > prob_E_I_connectivity) = 0;
% W_EI=W_EI/max(W_EI,[],'all');
W_EI = W_EI * InitialWeight; 

% W_IE(W_IE > prob_I_E_connectivity) = 0;
% W_IE=W_IE/max(W_IE,[],'all');
W_IE = W_IE * InitialWeight;


%% Running the network 
% Determining the speed of the animal running through the track
mean_dt = dt; sigma_dt = 0;              % time scale of updating activity
mean_v = (1); %%CHANGE DIRECTION 
sigma_v = 0;                    % mean and std speed of animal

dt_vec = normrnd(mean_dt, sigma_dt, 1, ceil(n_steps/mean_dt));
dt_vec(dt_vec < 0 ) = 0;
total_time_vec = cumsum(dt_vec);
v_vec0 = normrnd(mean_v, sigma_v, 1, ceil(n_steps/mean_dt)); 
v_vec = [0, v_vec0];
v_vec(v_vec < 0) = 0; %% CHANGE INEQUALITY

% Determining the path of the animal running through the track
positions = 0.5 * (v_vec(1:end-1) + v_vec(2:end)) .* dt_vec; % generating the positions at each time step from dt and the difference in velocities (trying to make it smooth)
positions = mod(cumsum(positions), n_pos);


pf_cell=pf_cell1;

RunSlowTm_STDP;

% Now you can save the results outside the parfor loop
    % Construct a unique filename for each iteration
    filename = sprintf('./SavedFiles/resultsFastSTDP1.mat');
    % Save the results
    save(filename, 'results');

        % Construct a unique filename for each iteration
    filename2 = './SavedFiles/W_InputEFastSTDP1.mat';
    % Save the results
    save(filename2, 'W_inputE');

    % Construct a unique filename for each iteration
    filename3 = './SavedFiles/PF_CellFastSTDP1.mat';
    % Save the results
    save(filename3, 'pf_cell1');

        % Construct a unique filename for each iteration
    filename4 = './SavedFiles/Mean_EIFastSTDP1.mat';
    % Save the results
    save(filename4, 'Mean_EI');

    % Construct a unique filename for each iteration
    filename5 = './SavedFiles/Mean_IEFastSTDP1.mat';
    % Save the results
    save(filename5, 'Mean_IE');

    % Construct a unique filename for each iteration
    filename6 = './SavedFiles/PosSTDP1.mat';
    % Save the results
    save(filename6, 'POS');

pf_cell=pf_cell2;

W_EE = 0;%rand(n_excit, n_excit); 
W_II = 0;%rand(n_inhib, n_inhib); 
W_EI = rand(n_excit, n_inhib); 
W_IE = rand(n_inhib, n_excit);

% W_EE(W_EE > prob_E_recurrent_connectivity) = 0;
% W_EE(diag(ones(1,n_excit)) == 1) = 0;       % no neurons can connect recurrently to themselves
% W_EE=W_EE/max(W_EE,[],'all');
% W_EE = W_EE * initial_weight_max_EE;
% 
% W_II(W_II > prob_I_I_connectivity) = 0;
% W_II(diag(ones(1,n_inhib)) == 1) = 0;       % no neurons can connect recurrently to themselves
% W_II=W_II/max(W_II,[],'all');
% W_II = W_II * initial_weight_max_II;


% W_EI(W_EI > prob_E_I_connectivity) = 0;
% W_EI=W_EI/max(W_EI,[],'all');
W_EI = W_EI * InitialWeight; 

% W_IE(W_IE > prob_I_E_connectivity) = 0;
% W_IE=W_IE/max(W_IE,[],'all');
W_IE = W_IE * InitialWeight;



RunSlowTm_STDP;

% Now you can save the results outside the parfor loop
    % Construct a unique filename for each iteration
    filename = sprintf('./SavedFiles/resultsFastSTDP2.mat');
    % Save the results
    save(filename, 'results');

        % Construct a unique filename for each iteration
    filename2 = './SavedFiles/W_InputEFastSTDP2.mat';
    % Save the results
    save(filename2, 'W_inputE');

    % Construct a unique filename for each iteration
    filename3 = './SavedFiles/PF_CellFastSTDP2.mat';
    % Save the results
    save(filename3, 'pf_cell2');

        % Construct a unique filename for each iteration
    filename4 = './SavedFiles/Mean_EIFastSTDP2.mat';
    % Save the results
    save(filename4, 'Mean_EI');

    % Construct a unique filename for each iteration
    filename5 = './SavedFiles/Mean_IEFastSTDP2.mat';
    % Save the results
    save(filename5, 'Mean_IE');

    % Construct a unique filename for each iteration
    filename6 = './SavedFiles/PosSTDP2.mat';
    % Save the results
    save(filename6, 'POS');

