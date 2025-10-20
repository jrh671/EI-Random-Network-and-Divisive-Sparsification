%% Initializing matrices & vectors

spike_mat_excitB = zeros(n_excit, length(total_time_vec)); 
spike_mat_inhibB = zeros(n_inhib, length(total_time_vec));
inhib_cumul_input = zeros(n_inhib, length(total_time_vec));
VmEB = zeros(n_excit, length(total_time_vec));
InputSpikesB = zeros(n_input, length(total_time_vec));
NetCurrentET = zeros(n_excit, length(total_time_vec));
NetCurrentIT = zeros(n_inhib, length(total_time_vec));
FrustrationET = zeros(n_excit, length(total_time_vec));
FrustrationIT = zeros(n_inhib, length(total_time_vec));
excit_spikesB = zeros(1, n_excit);                       
inhib_spikesB = zeros(1, n_inhib);                       
excit_spikes_2B = zeros(1, n_excit);                     
inhib_spikes_2B = zeros(1, n_inhib); 
excit_cum_inputB = zeros(1, n_excit);                    
inhib_cum_inputB = zeros(1, n_inhib); 
x_excit = zeros(1, n_excit);                            
x_inhib = zeros(1, n_inhib);
input_most_recent_fire_times_vecB = zeros(1, n_input) - 100; 
excit_most_recent_fire_times_vecB = zeros(1, n_excit) - 100;
inhib_most_recent_fire_times_vecB = zeros(1, n_inhib) - 100;
AveFrustIaT = zeros(1,length(total_time_vec));
AveFrustIbT = zeros(1,length(total_time_vec));
AveFrustET = zeros(1,length(total_time_vec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 W_EEf=zeros(n_excit,n_excit);
 W_EIf=ones(n_excit,n_inhib).*0.05;
 W_IEf=ones(n_inhib,n_excit).*0.05;
 
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
%     input_prob_firing_vec = input_prob_firing_vec;
    input_prob_firing_vec(input_prob_firing_vec < 0) = 0; 
    input_prob_firing_vec = (input_prob_firing_vec + alpha_P * pinknoise(n_input))* pf_rate * dt;
    input_prob_firing_vec(input_prob_firing_vec < 0) = 0; 

    
    % Finding which input neurons fired 
    coin_flips = rand(1, n_input);
    input_neurons_fired = find(coin_flips < input_prob_firing_vec & ...
        (total_time_vec(tt) - input_most_recent_fire_times_vecB) > input_refract_length);
    input_most_recent_fire_times_vecB(input_neurons_fired) = total_time_vec(tt);
    input_spikes = zeros(1, n_input); 
    input_spikes(input_neurons_fired) = 1; 
    InputSpikesB(:,tt)= input_spikes;
    % Updating excitatory cumulative input 
    excit_cum_inputB = excit_cum_inputB * exp(- dt) + input_spikes * W_inputE + ...
        excit_spikesB * W_EEf - inhib_spikesB * W_IEf + alpha_P * pinknoise(n_excit);     
    excit_cum_inputB((total_time_vec(tt) - excit_most_recent_fire_times_vecB) < excit_refract_length) = 0; % setting voltage equal to 0 for any neurons still in their absolute refractory period
    
    % Updating excitatory spikes 
    excit_most_recent_fire_times_vecB(excit_cum_inputB >= I_F_thresh) = total_time_vec(tt);
    excit_spikes_2B(excit_cum_inputB >= I_F_thresh) = 1;  % storing current timestep spikes
    excit_spikes_2B(excit_cum_inputB < I_F_thresh) = 0;
    excit_cum_inputB(excit_cum_inputB >= I_F_thresh) = 0; % setting the input back to 0 for the neurons that spiked        
    VmEB(:, tt) = excit_cum_inputB;
    spike_mat_excitB(:, tt) = excit_spikes_2B;
       
    % Updating inhibitory cumulative input 
    inhib_cum_inputB = inhib_cum_inputB * exp(-dt) ... 
    + excit_spikesB * W_EIf + alpha_P * pinknoise(n_inhib);     
    inhib_cum_inputB((total_time_vec(tt) - inhib_most_recent_fire_times_vecB) < inhib_refract_length) = 0;
    
    % Updating inhibitory spikes
    inhib_most_recent_fire_times_vecB(inhib_cum_inputB >= I_F_thresh) = total_time_vec(tt);
    inhib_spikes_2B(inhib_cum_inputB >= I_F_thresh) = 1;
    inhib_spikes_2B(inhib_cum_inputB < I_F_thresh) = 0;
    inhib_cum_inputB(inhib_cum_inputB >= I_F_thresh) = 0;     
    spike_mat_inhibB(:, tt) = inhib_spikes_2B;
    
    % Updating the spikes for the next times step
    inhib_spikesB = inhib_spikes_2B; 
    clear excit_spikesB
    excit_spikesB = excit_spikes_2B;
                            
end

FiringRateExb = SlidingAverage_s(spike_mat_excitB);

ComputePrimacySets;




