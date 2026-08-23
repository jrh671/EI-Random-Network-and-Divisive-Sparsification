if SlowTm==1      

UnitAdjInput = 1/mean_dt; 
UnitAdjInhibEI = 1/mean_dt; 
UnitAdjInhibIE = 1/mean_dt;

C=1;
gL=1;
Tm=(C/gL); %Resulting units are in seconds
VL=0;
VLeakE = VL*ones(1,n_excit); %Leak Potential
VLeakI = VL*ones(1,n_inhib); %Leak Potential
VL = 0;
VrE = 0;
VrI = 0;

%Integrate and Fire Parameters
        I_F_threshE = 1; 
        I_F_threshI = 1; 

        input_refract_length = 0;
        excit_refract_length = 0.02;
        inhib_refract_length = 0.02;

        %Connection Scaling Factors
        eta = 0.0005;
        prob_pf_center = 0.20; 
        prob_E_recurrent_connectivity = 0.25;
        prob_I_E_connectivity = 0.30;
        prob_E_I_connectivity = 0.175;
        prob_I_I_connectivity = 0.50;

        %Connection Limits
        initial_weight_max_input = 0.05;
        initial_weight_max_EE = 0.05 / prob_E_recurrent_connectivity;
        initial_weight_max_EI = 0.05 / prob_I_E_connectivity;
        initial_weight_max_IE = 0.05 / prob_E_I_connectivity;
        initial_weight_max_II = 0.05 / prob_I_I_connectivity;
        W_upper_limit = 0.5;
        dt=mean_dt;

        pf_rate = 12; % Input neuron rate (in Hz) at center of pf
        pf_width_ = .3162; % Controls Gaussian pf width pf_width = (1/Sigma^2)

else

    %% Parameters 

MaxWeight=0.5;
%% Adjusting Input, IE and EI resistance to current
UnitAdjInput = .00015*10; 
UnitAdjInhibEI = .015*15; 
UnitAdjInhibIE = .015/10;

prob_pf_center = 0.20; 
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
pf_rate = 25;                            % input neuron rate (in Hz) at center of pf
pf_width_ = 2.4495;                              % controls the Gaussian pf width (note it's actually the inverse) 
SynConst = 1; %Synaptic Time Constant for Plasticity
corr_window_size = 0.25;                    % size of window used for smoothing activity for correlation analysis
dt=mean_dt;
Tsyn=1;

EI_plast=1;
II_plast=0;
IE_plast=1;
EE_plast=0;
eta=0.0005;


end

        alpha_P = 0; %Noise Input
        alpha_Pe = 0; %Noise E
        alpha_Pi = 0; %Noise I

        %% Initialize weights and place fields
        InputStrength = 0.05;
        W_inputE = InputStrength * rand(n_input, n_excit);
        W_EE = 0;
        W_EI = 0.05 * ones(n_excit, n_inhib);
        W_IE = 0.05 * ones(n_inhib, n_excit);
        W_II = 0;
