addpath('./Helper_Functions');
addpath('./SavedFiles');

%% Adjust Save Directory Below (Lines 280-308 and Lines 493-517)

rng(0)

%% Variables to predefine outputs (If changed, make sure to also changes within parfor!)
NumIterations=1; 
NumRepetitions=1; 
NumNeurons=500;
NumInputs=10000;
NumPos = 30;
NumLaps=5;
ProbPF = 0.2;
InitialWeight=0.005;

pf_cell = GET_TRACK_PFs(NumPos,ProbPF, NumInputs);

      
n_laps = 5; 
n_pos = 30;

MaxWeight=0.5;

%% Adjusting Input, IE and EI resistance to current
UnitAdjInput = .00015; 
UnitAdjInhibEI = .075; 
UnitAdjInhibIE = .075/20;
 
%% Parameters 
n_input = 10000;                             % number of tuned input neurons                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
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
gL=.000025; %Conductance u(mS) where mS = (C/mV)   | Note S.I. S = (C/V)

Tm=(C/gL); %Resulting units are in seconds


n_steps = n_laps * n_pos;                   % number of steps on the track 
input_refract_length = 0;                   % input neuron refractory length 
excit_refract_length = 0.002;                % excit refractory length
inhib_refract_length = 0.002;                % inhibitory refractory length
W_upper_limit = MaxWeight;                        % maximum weight we allow the E to I weights to take
alpha_P = 0;                             % pink noise strength
pf_rate = 25;                                % input neuron rate (in Hz) at center of pf
pf_width = 1/6;                              % controls the Gaussian pf width (note it's actually the inverse) 
SynConst = 1; %Synaptic Time Constant for Plasticity
corr_window_size = 0.25;                    % size of window used for smoothing activity for correlation analysis
dt=0.0005;
Tsyn=.020;

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


%% Initializing matrices & vectors to store various information about the
% dynamics

spike_mat_excit = zeros(n_excit, length(total_time_vec)); % stores excitatory spikes 
spike_mat_inhib = zeros(n_inhib, length(total_time_vec));
VmI = zeros(n_inhib, length(total_time_vec));
VmE = zeros(n_excit, length(total_time_vec));
Mean_EI = zeros(length(total_time_vec),1);
Mean_IE = zeros(length(total_time_vec),1);
NetCurrentET = zeros(n_excit, length(total_time_vec));
NetCurrentIT = zeros(n_inhib, length(total_time_vec));
FrustrationET = zeros(n_excit, length(total_time_vec));
FrustrationIT = zeros(n_inhib, length(total_time_vec));
excit_spikes = zeros(1, n_excit);                       % stores the spiking at each time step of the excitatory neurons 
inhib_spikes = zeros(1, n_inhib);                       
excit_spikes_2 = zeros(1, n_excit);                     % buffer that stores the t+1 activity so we can still keep the t activity
inhib_spikes_2 = zeros(1, n_inhib); 
excit_cum_input = zeros(1, n_excit);                    % stores the voltage of the excitatory neurons
inhib_cum_input = zeros(1, n_inhib); 
x_excit = zeros(1, n_excit);                            % sotres the exponentials used in calculating the STDP changes
x_inhib = zeros(1, n_inhib);
input_most_recent_fire_times_vec = zeros(1, n_input) - 100; % stores the most recent firing times of the input neurons 
excit_most_recent_fire_times_vec = zeros(1, n_excit) - 100;
inhib_most_recent_fire_times_vec = zeros(1, n_inhib) - 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Program Stimulation
% Running the network on the linear track 
for tt = 1:length(total_time_vec)

    if mod(tt,3/mean_dt)==1
        ceil(tt/(3/mean_dt))
    end



    % Getting the position and time step
    pos = positions(tt);
    dt = dt_vec(tt);
    
    % Finding how far away the animal is from all of the rf centers and
    % calculating the probability of each of the input neurons firing
    input_prob_firing_vec = zeros(1, n_input); 
    for ii = 1:n_input
        rfs = pf_cell{ii}; 
        possible_distance_vec = [abs(pos - rfs); (n_pos) - abs(pos - rfs)];
        distance_vec = min(possible_distance_vec); 
        input_prob_firing_vec(ii) = sum(exp(-distance_vec.^2 * pf_width)); 
    end
    input_prob_firing_vec = input_prob_firing_vec + alpha_P * pinknoise(n_input);
    input_prob_firing_vec(input_prob_firing_vec < 0) = 0; 
    input_prob_firing_vec = input_prob_firing_vec * pf_rate * dt; 
    
    % Finding which input neurons fired 
    coin_flips = rand(1, n_input);
    input_neurons_fired = find(coin_flips < input_prob_firing_vec & ...
        (total_time_vec(tt) - input_most_recent_fire_times_vec) > input_refract_length);
    input_most_recent_fire_times_vec(input_neurons_fired) = total_time_vec(tt);
    input_spikes = zeros(1, n_input); 
    input_spikes(input_neurons_fired) = 1; 
    
  
  

% Updating excitatory cumulative input 
    excit_cum_input = excit_cum_input + ((1)/C)*input_spikes * W_inputE*UnitAdjInput*((dt)) + ...
       (1/C)*excit_spikes * W_EE*dt - ((gL/C))*(excit_cum_input-VLeakE)*dt - (1/C)*inhib_spikes * W_IE*UnitAdjInhibIE*(dt) + alpha_P * pinknoise(n_excit);     
    excit_cum_input((total_time_vec(tt) - excit_most_recent_fire_times_vec) < excit_refract_length) = VrE; % setting voltage equal to 0 for any neurons still in their absolute refractory period
    
    % Updating excitatory spikes 
    excit_most_recent_fire_times_vec(excit_cum_input >= I_F_threshE) = total_time_vec(tt);
    excit_spikes_2(excit_cum_input >= I_F_threshE) = 1;  % storing current timestep spikes
    excit_spikes_2(excit_cum_input < I_F_threshE) = 0;
    excit_cum_input(excit_cum_input >= I_F_threshE) = VrE; % setting the input back to 0 for the neurons that spiked        
    spike_mat_excit(:, tt) = excit_spikes_2;
       
    % Updating inhibitory cumulative input 
    inhib_cum_input = inhib_cum_input - (gL/C)*(inhib_cum_input-VLeakI)*dt + ...
        (1/C)*excit_spikes * W_EI*UnitAdjInhibEI*(dt) - (1/C)*inhib_spikes * W_II*dt + alpha_P * pinknoise(n_inhib);     
    inhib_cum_input((total_time_vec(tt) - inhib_most_recent_fire_times_vec) < inhib_refract_length) = VrI;
    
    % Updating inhibitory spikes
    inhib_most_recent_fire_times_vec(inhib_cum_input >= I_F_threshI) = total_time_vec(tt);
    inhib_spikes_2(inhib_cum_input >= I_F_threshI) = 1;
    inhib_spikes_2(inhib_cum_input < I_F_threshI) = 0;
    inhib_cum_input(inhib_cum_input >= I_F_threshI) = VrI;     

    % Updating the spikes for the next times step
    inhib_spikes = inhib_spikes_2; 
    excit_spikes = excit_spikes_2;

    % Updating the inhibitory/excitatory exponentials for STDP
    x_excit = exp(-dt/Tsyn) * x_excit;
    x_excit(excit_spikes == 1) = 1 + x_excit(excit_spikes == 1); 
    x_inhib = exp(-dt/Tsyn) * x_inhib;
    x_inhib(inhib_spikes == 1) = 1 + x_inhib(inhib_spikes == 1);
             
    
    % Updating weights of the network. Remember that the input weights are 
    % not plastic but the other weights are 
    d_W_EI = (inhib_spikes' * x_excit)' - excit_spikes' * x_inhib;
    d_W_EI(d_W_EI > 1) = 1; d_W_EI(d_W_EI < -1) = -1; % making the max change in weight be eta
    d_W_IE = (excit_spikes' * x_inhib)' + inhib_spikes' * x_excit - 0.0001 * ones(size(W_IE)); % shifting the curve by .1%
    d_W_IE(d_W_IE > 1) = 1; d_W_IE(d_W_IE < -1) = -1;
    d_W_EE = x_excit' * excit_spikes - excit_spikes' * x_excit; 
    d_W_EE(d_W_EE > 1) = 1; d_W_EE(d_W_EE < -1) = -1;
    d_W_II = x_inhib' * inhib_spikes - inhib_spikes' * x_inhib; 
    d_W_II(d_W_II > 1) = 1; d_W_II(d_W_II < -1) = -1;
     

    W_EI = W_EI + EI_plast * eta * d_W_EI;    
    W_IE = W_IE + IE_plast * eta * d_W_IE; 
    W_EE = W_EE + EE_plast * eta * d_W_EE; 
    W_II = W_II + II_plast * eta * d_W_II; 

    % Enforcing weight boundaries 
    W_EE(W_EE < 0 ) = 0; 
    W_EI(W_EI < 0) = 0; 
    W_IE(W_IE < 0) = 0;
    W_II(W_II < 0) = 0;
    W_EI(W_EI > W_upper_limit) = W_upper_limit; 
    W_IE(W_IE > W_upper_limit) = W_upper_limit; 
    W_EE(W_EE > W_upper_limit) = W_upper_limit; 
    W_II(W_II > W_upper_limit) = W_upper_limit; 
     
    Mean_EI(tt) = mean(W_EI,'all');
    Mean_IE(tt) = mean(W_IE,'all');
    
end


Rates=SlidingAverage_s(spike_mat_excit,mean_dt/mean_v);
results= Rates;

threshold=0.975;

 
% Now you can save the results outside the parfor loop
    % Construct a unique filename for each iteration
    filename = sprintf('./Directory_To/resultsFastSTDP1.mat');
    % Save the results
    save(filename, 'results');

        % Construct a unique filename for each iteration
    filename2 = './Directory_To/W_InputEFastSTDP1.mat';
    % Save the results
    save(filename2, 'W_inputE');

    % Construct a unique filename for each iteration
    filename3 = './Directory_To/PF_CellFastSTDP1.mat';
    % Save the results
    save(filename3, 'pf_cell');

        % Construct a unique filename for each iteration
    filename4 = './Directory_To/Mean_EIFastSTDP1.mat';
    % Save the results
    save(filename4, 'Mean_EI');

    % Construct a unique filename for each iteration
    filename5 = './Directory_To/Mean_IEFastSTDP1.mat';
    % Save the results
    save(filename5, 'Mean_IE');


