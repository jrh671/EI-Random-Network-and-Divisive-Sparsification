function [spike_mat_excit,positions,DownstreamNeuronOutputs,DownstreamNeuronOutput,ReverseData,FreezeData, FreezeNowRecord,DirectionChangeRecord,FR] = RRN_SensorimotorRecallNetwork(n_input,n_excit,n_inhib,W_inputE,pf_cell1,pf_cell2,GoalNeurons,PrimacyNeurons)

%% Contextual Switch Times (Note Analysis Calculations assume even 1/3rds for A,B,A contexts.
PFSwitchOn = 301; %Begin Context 2 (randomly reassign input filters)
PFSwitchOff = 600; %Return to Context 1 (reset to initial input filters)

% Recurrent weights
W_EEf=zeros(n_excit,n_excit);
W_EIf=ones(n_excit,n_inhib).*0.05;
W_IEf=ones(n_inhib,n_excit).*0.05;

%% Track Setip
n_laps = 30; %Each Lap = n_pos sec
n_pos = 30; 

%% Parameters 
n_steps = n_laps * n_pos;                   % number of steps on the track 
I_F_thresh = 1;                           % integrate and fire threshold 

input_refract_length = 0;                   % input neuron refractory length 
excit_refract_length = 0.02;                % excit refractory length
inhib_refract_length = 0.02;                % inhibitory refractory length
alpha_P = 0;                             % pink noise strength
pf_rate = 8;                                % input neuron rate (in Hz) at center of pf
pf_width = 10;                              % controls the Gaussian pf width (note it's actually the inverse) 10

% Determining the speed of the animal running through the track
mean_dt = 0.005; sigma_dt = 0;              % time scale of updating activity
mean_v = 1; 
sigma_v = 0;       

dt_vec = normrnd(mean_dt, sigma_dt, 1, ceil(n_steps/mean_dt));
dt_vec(dt_vec < 0 ) = 0;
total_time_vec = cumsum(dt_vec);
v_vec0 = normrnd(mean_v, sigma_v, 1, ceil(n_steps/mean_dt)); 
v_vec = [0, v_vec0];
v_vec(v_vec < 0) = 0; 


%% Initializing matrices & vectors to store data
spike_mat_excit = zeros(n_excit, length(total_time_vec)); % stores excitatory spikes 
spike_mat_inhib = zeros(n_inhib, length(total_time_vec));
excit_spikes = zeros(1, n_excit);                       % stores the spiking at each time step of the excitatory neurons 
inhib_spikes = zeros(1, n_inhib);                       
excit_spikes_2 = zeros(1, n_excit);                     % buffer that stores the t+1 activity so we can still keep the t activity
inhib_spikes_2 = zeros(1, n_inhib); 
excit_cum_input = zeros(1, n_excit);                    % stores the voltage of the excitatory neurons
inhib_cum_input = zeros(1, n_inhib); 
input_most_recent_fire_times_vec = zeros(1, n_input) - 100; % stores the most recent firing times of the input neurons 
excit_most_recent_fire_times_vec = zeros(1, n_excit) - 100;
inhib_most_recent_fire_times_vec = zeros(1, n_inhib) - 100;
Input_Spikes = zeros(n_input,length(total_time_vec));
DownstreamNeuronOutput = zeros(size(spike_mat_excit(1,:)));
DownstreamNeuronOutputs = zeros(25,length(spike_mat_excit(1,:)));

% Time Windows
TotalTimeBins = length(total_time_vec);
freeze_duration = 5 / mean_dt;  % Duration of freeze in time bins
decision_interval = ceil(1 / mean_dt);  % Time bins between decision computations
FreezeNowRecord = zeros(1, TotalTimeBins);
FreezeNow = 0;
DirectionChangeRecord = zeros(1, TotalTimeBins);
DirectionChangeNow = 0;
positions = zeros(1, TotalTimeBins);
last_position_before_freeze = 0;
freeze_on_this_lap = false;
direction_change_on_this_lap = false;
time_to_freeze = 0;
time_per_lap = n_pos / mean_v / mean_dt;  % Time bins per lap
lap_start = 1;
decision_time = 5*decision_interval;  % Initialize to the first decision time
last_three_seconds_activity = zeros(500,ceil(3/mean_dt));
% Initialize last_reversal_time at the start of your code
last_reversal_time = -Inf; % Set it to -Inf initially to allow the first reversal
direction = 1; %Initialize Direction

%% Begin Recall Simulation
 
for tt = 1:ceil(TotalTimeBins)
 
  %% Sliding windows to check for downstream activations by tagged neurons
    if tt == decision_time

        % Determine downstream if Primacy Sets were activated
        actions = RRN_DownstreamNeurons(last_three_seconds_activity,PrimacyNeurons); 
        % Determine downstream if Goal Cells were activated
        direction_decision = RRN_DownstreamNeuron(last_three_seconds_activity, GoalNeurons); 

        % Store downstream output 
        DownstreamNeuronOutputs(:,tt) = actions;
        DownstreamNeuronOutput(tt) = direction_decision;
        
        action = sum(actions);
        direction_change = direction_decision; % Current direction of motion

        % Freeze if Primacy reactivation occurs
        if action > 1 && time_to_freeze <= 0 
            'Freeze';
            time_to_freeze = freeze_duration;  
            FreezeNow = FreezeNow+1;
            FreezeNowRecord(tt) = FreezeNow;
            freeze_on_this_lap = true;
        end
       
        % % Change direction if Goal reactivation occurs
        % if direction_change == 1 && direction_change_on_this_lap == false
        %     'Reversal';
        %     direction = -direction;  % flip direction
        %     direction_change_on_this_lap = true;
        %     DirectionChangeNow = DirectionChangeNow+1;
        %     DirectionChangeRecord(tt) = DirectionChangeNow;
        % 
        % end


        % Change direction if Goal reactivation occurs and meets the 5-second condition
        if direction_change == 1 && direction_change_on_this_lap == false
            % Check if enough time has passed since the last reversal
            if tt - last_reversal_time >= (1/mean_dt) * 10
                'Reversal';
                direction = -direction;  % flip direction
                direction_change_on_this_lap = true;
                DirectionChangeNow = DirectionChangeNow + 1;
                DirectionChangeRecord(tt) = DirectionChangeNow;
        
                % Update the last reversal time
                last_reversal_time = tt;
            else
                'Reversal blocked due to time condition'; % Optional debug message
            end
        end


        decision_time = decision_time + decision_interval;  % Schedule the next decision time
    end

    %% Update position and action sequences according to downstream output
    if tt>1

        if time_to_freeze > 0
    
            % If the action was to freeze, the position stays the same
            positions(tt) = positions(tt - 1);
            time_to_freeze = time_to_freeze - 1;
    
        elseif time_to_freeze <= 0 
    
            % Otherwise, the position updates normally according to the speed and direction
            positions(tt) = positions(tt - 1) + direction * mean_v * dt_vec(tt);
            positions(tt) = mod(positions(tt), n_pos);  % To keep the position within the track limits
        end

        % Direction changes respecting circular boundary conditions
        if mod(positions(tt), n_pos) < mod(positions(tt-1), n_pos) && direction == 1 && direction_change_on_this_lap == true
            freeze_on_this_lap = false;
            direction_change_on_this_lap = false;
            lap_start = tt;
    
            elseif mod(positions(tt), n_pos) > mod(positions(tt-1), n_pos) && direction == -1 && direction_change_on_this_lap == true
                freeze_on_this_lap = false;
                direction_change_on_this_lap = false;
                lap_start = tt;
        end

    else %Second time point and beyond
        positions(tt) = direction * mean_v * dt_vec(tt);
        positions(tt) = mod(positions(tt), n_pos);  % To keep the position within the track limits
    end

    %% Simulate Positional Input Activity According To "Online" Actions
    pos = positions(tt);
    dt = dt_vec(tt);
    input_prob_firing_vec = zeros(1, n_input); 

    % Finding distance from all place fields for a given input
    for ii = 1:n_input
        if tt>PFSwitchOn*200 && tt<PFSwitchOff*200
            clear rfs
                rfs = pf_cell1{ii}; 
        else
            clear rfs
                rfs = pf_cell2{ii}; 
            
        end
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
    Input_Spikes(:,tt) = input_spikes;
    
    %% EI Network

    % Updating excitatory cumulative input 
    excit_cum_input = excit_cum_input * exp(- dt) + input_spikes * W_inputE + ...
        excit_spikes * W_EEf - inhib_spikes * W_IEf + alpha_P * pinknoise(n_excit);     
    excit_cum_input((total_time_vec(tt) - excit_most_recent_fire_times_vec) < excit_refract_length) = 0; % setting voltage equal to 0 for any neurons still in their absolute refractory period
    
    % Updating excitatory spikes 
    excit_most_recent_fire_times_vec(excit_cum_input >= I_F_thresh) = total_time_vec(tt);
    excit_spikes_2(excit_cum_input >= I_F_thresh) = 1;  % storing current timestep spikes
    excit_spikes_2(excit_cum_input < I_F_thresh) = 0;
    excit_cum_input(excit_cum_input >= I_F_thresh) = 0; % setting the input back to 0 for the neurons that spiked        
    VmE(:, tt) = excit_cum_input;
    spike_mat_excit(:, tt) = excit_spikes_2;
      
        % Updating inhibitory cumulative input 
    inhib_cum_input = inhib_cum_input * exp(-dt) ...
        + excit_spikes * W_EIf + alpha_P * pinknoise(n_inhib);     
    inhib_cum_input((total_time_vec(tt) - inhib_most_recent_fire_times_vec) < inhib_refract_length) = 0;
    
    % Updating inhibitory spikes
    inhib_most_recent_fire_times_vec(inhib_cum_input >= I_F_thresh) = total_time_vec(tt);
    inhib_spikes_2(inhib_cum_input >= I_F_thresh) = 1;
    inhib_spikes_2(inhib_cum_input < I_F_thresh) = 0;
    inhib_cum_input(inhib_cum_input >= I_F_thresh) = 0;     
    VmI(:,tt) = inhib_cum_input;
    spike_mat_inhib(:, tt) = inhib_spikes_2;

    % Updating the spikes for the next times step
    inhib_spikes = inhib_spikes_2; 
    excit_spikes = excit_spikes_2;
    
    % Record the spiking activity
    spiking_activity(:,tt) = spike_mat_excit(:,tt); % Update this with your own code
    
    % Update the last 3 seconds activity
    last_three_seconds_activity = [last_three_seconds_activity(:,2:end), spiking_activity(:,tt)];
  
    %% End EI Network
                     
end

%% Save Results
FR = SlidingAverage_s(spike_mat_excit);
FRi = SlidingAverage_s(spike_mat_inhib);

% Number of bins per epoch
EpochBins = TotalTimeBins / 3;

% Collect cumulative freeze results into three separate epochs
Freeze1 = max(FreezeNowRecord(1:EpochBins));  % First epoch value
Freeze2 = max(FreezeNowRecord(1:2*EpochBins)) - max(FreezeNowRecord(1:EpochBins));  % Second epoch value
Freeze3 = max(FreezeNowRecord(1:TotalTimeBins)) - max(FreezeNowRecord(EpochBins+1:2*EpochBins));  % Third epoch value

% Collect cumulative reversal results into three separate epochs
Reverse1 = max(DirectionChangeRecord(1:EpochBins));  % First epoch value
Reverse2 = max(DirectionChangeRecord(1:2*EpochBins)) - max(DirectionChangeRecord(1:EpochBins));  % Second epoch value
Reverse3 = max(DirectionChangeRecord(1:TotalTimeBins)) - max(DirectionChangeRecord(EpochBins+1:2*EpochBins));  % Third epoch value

% Store results
FreezeData = [Freeze1, Freeze2, Freeze3];
ReverseData = [Reverse1, Reverse2, Reverse3];

ReverseData
% FreezeData(FreezeData < 0) = 0;


