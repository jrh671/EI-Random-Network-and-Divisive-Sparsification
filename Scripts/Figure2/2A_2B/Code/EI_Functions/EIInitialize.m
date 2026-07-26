
%% Track Setip
n_laps = 30; %Each Lap = n_pos sec
n_pos = 30;
InputStrength = 1;

if J==1
CollectPositions = zeros(n_laps*n_pos,numSims);
end

%% WeightI
EI = 0.05; %EI Weight .05
IE = 0.05; %IE Weight .05

%% Weight EE
W_EE = 0;%Or Point 2

On = 1; 
Off = 0;



%% Goal Cell Functions
rng(J+1e2) %Same seed as Context 2
special_position = randperm(n_pos,1);
SpecialPositions(J) = special_position;
Fr= 2.75;
prop_special_neurons = 0.2;


%% Disinhibition Grouping
GroupDegree = 1;

n_input = 1000;                             % number of tuned input neurons                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
n_excit = 500;                              % number of excitatory neurons

if GroupDegree == 1
n_inhib = 50;                               % number of inhibititory neurons 
end

if GroupDegree == 2
n_inhib = 100;                               % number of inhibititory neurons 
end

if GroupDegree == 3
n_inhib = 150;                               % number of inhibititory neurons 
end

if GroupDegree == 4
n_inhib = 200;                               % number of inhibititory neurons 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NovelPF = On; 
NovelInput = On;


%% Parameters 
n_steps = n_laps * n_pos;                   % number of steps on the track 
eta = 0.0005;  
prob_pf_center = 0.20; 
prob_E_recurrent_connectivity = 0.25;       % probability of any given excitatory neuron will connect to another excitatory neuron
prob_I_E_connectivity = 0.30;               % probability of any inhibitory neuron connecting to an excitatory neuron
prob_E_I_connectivity = 0.175;
prob_I_I_connectivity = .50;              % probability of any I neuron connecting to an inhibitory neuron                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ; %probability of any excitatory neuron connecting to an inhibitory neuron
I_F_thresh = 1;                           % integrate and fire threshold 
test_time = 10;                             % number of laps used to compute metrics

input_refract_length = 0;                   % input neuron refractory length 
excit_refract_length = 0.02;                % excit refractory length
inhib_refract_length = 0.02;                % inhibitory refractory length
initial_weight_max_input = 0.05;            % initial weight maximum input to E
initial_weight_max_EE = 0.05 / prob_E_recurrent_connectivity; % initial weight maximum for E to E
initial_weight_max_EI = 0.05 / prob_I_E_connectivity;  % initial weight maximum for E to I
initial_weight_max_IE = 0.05 / prob_E_I_connectivity; % initial weight maximum for I to E
initial_weight_max_II = 0.05/ prob_I_I_connectivity;
W_upper_limit = 0.5;                        % maximum weight we allow the E to I weights to take
alpha_P = 0;                             % pink noise strength
pf_rate = 8;                                % input neuron rate (in Hz) at center of pf
pf_width = 10;                              % controls the Gaussian pf width (note it's actually the inverse) 10
corr_window_size = 0.25;                    % size of window used for smoothing activity for correlation analysis

%% Initiliazing the input neurons 
            

% Note: W_inputE is not plastic
if NovelInput == 1
    rng(0) %Ensure the same network is always used
W_inputE = InputStrength*rand(n_input, n_excit) * 0.05;

end



 W_EEf=zeros(n_excit,n_excit);
 W_IIf=zeros(n_inhib,n_inhib);
 W_EIf=ones(n_excit,n_inhib).*EI;
 W_IEf=ones(n_inhib,n_excit).*IE;
