%% Variables to predefine outputs (If changed, make sure to also changes within parfor!)
NumIterations=21; 
NumRepetitions=10; 
NumNeurons=500;
NumInputs=1000;
NumPos = 30;
NumLaps=5;
ProbPF = 0.2;
ScaleInputs=0.05;

results = cell(1 * NumIterations ,1);

pf_cellC = GET_TRACK_PFs(NumPos,ProbPF, NumInputs);
W_inputECell = cell(1,1);
W_inputECell{1} = rand(NumInputs, NumNeurons) * ScaleInputs;

for t=1:NumIterations*NumRepetitions
      results{t}=zeros(NumNeurons,NumPos*NumLaps);
end


parpool(32);

parfor idx = 1:(NumIterations*NumRepetitions)
% for idx = 10*21
    % Calculate P and N based on idxd
    sample = mod(idx-1, 10) + 1;
    Ohm = ceil(idx / 10);
%     Instance= ceil(idx / 210);
    W_inputElocal=W_inputECell;
    W_inputE=W_inputElocal{1};
    pf_cell=pf_cellC;


n_laps = 5; 
n_pos = 30;

%% Input, EI and IE Resistance to Current
InputWeight=1; 
UnitAdjInput = 1; 
UnitAdjInhibEI = 1; 
UnitAdjInhibIE = 1;

n_laps = n_laps+1; % Voltage Equilibrium Initial Lap.
 
%% Parameters 
n_input = 1000;                             % number of tuned input neurons                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
n_excit = 500;                              % number of excitatory neurons
n_inhib = 50;                               % number of inhibititory neurons 

%% Parameters 

prob_E_recurrent_connectivity = 0.25;       % probability of any given excitatory neuron will connect to another excitatory neuron
prob_I_E_connectivity = 0.30;               % probability of any inhibitory neuron connecting to an excitatory neuron
prob_E_I_connectivity = 0.175;
prob_I_I_connectivity = .50;              % probability of any I neuron connecting to an inhibitory neuron                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ; %probability of any excitatory neuron connecting to an inhibitory neuron
I_F_threshE = 1;                           % integrate and fire threshold 
I_F_threshI = 1;                           % integrate and fire threshold 

VL = 0;
VrE = 0;
VrI = 0;
VLeakE = VL*ones(1,n_excit); %Leak Potential
VLeakI = VL*ones(1,n_inhib); %Leak Potential

C=1; %Arbitrary units
gL=1; %Arbitrary units

Tm=(C/gL); %Seconds Units


n_steps = n_laps * n_pos;                   % number of steps on the track 
input_refract_length = 0;                   % input neuron refractory length 
excit_refract_length = 0.02;                % excit refractory length
inhib_refract_length = 0.02;                % inhibitory refractory length
alpha_P = 0;                             % pink noise strength
pf_rate = 12;                                % input neuron rate (in Hz) at center of pf
pf_width = 10;                              % controls the Gaussian pf width (note it's actually the inverse) 
dt=0.005;

%% Initiliazing the input neurons 
% Setting up the weight matrices

W_inputE= InputWeight*W_inputE;

%% II Weights %%
W_EE = 0*ones(n_excit,n_excit);

W_EI = ((Ohm-1)*0.015)*ones(n_excit,n_inhib);

W_IE = ((Ohm-1)*0.015)*ones(n_inhib,n_excit);

W_II = 0*rand(n_inhib,n_inhib);

W_II(1:n_inhib/2,1:n_inhib/2) = 0;
W_II(1+n_inhib/2:end,1+n_inhib/2:end) = 0;



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
    excit_cum_input = excit_cum_input * exp(-dt) + ((1)/C)*input_spikes * W_inputE*UnitAdjInput + ...
       (1/C)*excit_spikes * W_EE - (1/C)*inhib_spikes * W_IE*UnitAdjInhibIE + alpha_P * pinknoise(n_excit);     
    excit_cum_input((total_time_vec(tt) - excit_most_recent_fire_times_vec) < excit_refract_length) = VrE; % setting voltage equal to 0 for any neurons still in their absolute refractory period
    
    % Updating excitatory spikes 
    excit_most_recent_fire_times_vec(excit_cum_input >= I_F_threshE) = total_time_vec(tt);
    excit_spikes_2(excit_cum_input >= I_F_threshE) = 1;  % storing current timestep spikes
    excit_spikes_2(excit_cum_input < I_F_threshE) = 0;
    excit_cum_input(excit_cum_input >= I_F_threshE) = VrE; % setting the input back to 0 for the neurons that spiked        
    spike_mat_excit(:, tt) = excit_spikes_2;
 
       
       
    % Updating inhibitory cumulative input 
    inhib_cum_input = inhib_cum_input * exp(-dt) + ...
        (1/C)*excit_spikes * W_EI*UnitAdjInhibEI - (1/C)*inhib_spikes * W_II + alpha_P * pinknoise(n_inhib);     
    inhib_cum_input((total_time_vec(tt) - inhib_most_recent_fire_times_vec) < inhib_refract_length) = VrI;
    
    % Updating inhibitory spikes
    inhib_most_recent_fire_times_vec(inhib_cum_input >= I_F_threshI) = total_time_vec(tt);
    inhib_spikes_2(inhib_cum_input >= I_F_threshI) = 1;
    inhib_spikes_2(inhib_cum_input < I_F_threshI) = 0;
    inhib_cum_input(inhib_cum_input >= I_F_threshI) = VrI;     


    % Updating the spikes for the next times step
    inhib_spikes = inhib_spikes_2; 
    excit_spikes = excit_spikes_2;
       

                     
end
 

Rates=SlidingAverage_s(spike_mat_excit,mean_dt/mean_v);
results{idx}= Rates(:,31:end);

end


% Now you can save the results outside the parfor loop
    % Construct a unique filename for each iteration
    filename = sprintf('/mnt/home/jhurtado/ceph/resultsSlow.mat');
    % Save the results
    save(filename, 'results');

        % Construct a unique filename for each iteration
    filename2 = '/mnt/home/jhurtado/ceph/W_InputESlow.mat';
    % Save the results
    save(filename2, 'W_inputECell');

    % Construct a unique filename for each iteration
    filename3 = '/mnt/home/jhurtado/ceph/PF_CellSlow.mat';
    % Save the results
    save(filename3, 'pf_cellC');

% Close the parallel pool
delete(gcp('nocreate'));
