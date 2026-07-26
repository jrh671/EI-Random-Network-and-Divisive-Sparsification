for tt = 1:length(total_time_vec)

    pos = positions(tt);
    dt = dt_vec(tt);
    
    % Finding how far away the animal is from all of the rf centers and calculating the probability of each of the input neurons firing
    input_prob_firing_vec = zeros(1, n_input); 
    
    for ii = 1:n_input
            rfs = pf_cell{ii};
    
        possible_distance_vec = [abs(pos - rfs); (n_pos) - abs(pos - rfs)];
        distance_vec = min(possible_distance_vec); 
        input_prob_firing_vec(ii) = sum(exp(-(distance_vec.^2) * pf_width));
    end
    
   if tt > ceil(n_pos/mean_dt)
        % Modify firing rates at the special position
        if ceil(pos) == ceil(special_position)
            for ii = 1:n_input
                if ismember(ceil(special_position), ceil(pf_cell{ii}))
                    input_prob_firing_vec(ii) = input_prob_firing_vec(ii) * Fr;
                end
            end
        end
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
    Input_Spikes(tt,:) = input_spikes;
    
    % Updating excitatory cumulative input 
    excit_cum_inputf = excit_cum_inputf * exp(- dt) + input_spikes * W_inputE + ...
        excit_spikesf * W_EEf - inhib_spikesf * W_IEf + alpha_P * pinknoise(n_excit);     
    excit_cum_inputf((total_time_vec(tt) - excit_most_recent_fire_times_vecf) < excit_refract_length) = 0; % setting voltage equal to 0 for any neurons still in their absolute refractory period
    
    % Updating excitatory spikes 
    excit_most_recent_fire_times_vecf(excit_cum_inputf >= I_F_thresh) = total_time_vec(tt);
    excit_spikes_2f(excit_cum_inputf >= I_F_thresh) = 1;  % storing current timestep spikes
    excit_spikes_2f(excit_cum_inputf < I_F_thresh) = 0;
    excit_cum_inputf(excit_cum_inputf >= I_F_thresh) = 0; % setting the input back to 0 for the neurons that spiked        
    excit_cumul_inputf(:, tt) = excit_cum_inputf;
    spike_mat_excitf(:, tt) = excit_spikes_2f;
       
    % Updating inhibitory cumulative input 
    inhib_cum_inputf = inhib_cum_inputf * exp(-dt) - inhib_spikesf*W_IIf ...
        + excit_spikesf * W_EIf + alpha_P * pinknoise(n_inhib);     
    inhib_cum_inputf((total_time_vec(tt) - inhib_most_recent_fire_times_vecf) < inhib_refract_length) = 0;
    
    % Updating inhibitory spikes
    inhib_most_recent_fire_times_vecf(inhib_cum_inputf >= I_F_thresh) = total_time_vec(tt);
    inhib_spikes_2f(inhib_cum_inputf >= I_F_thresh) = 1;
    inhib_spikes_2f(inhib_cum_inputf < I_F_thresh) = 0;
    inhib_cum_inputf(inhib_cum_inputf >= I_F_thresh) = 0;     
    inhib_cumul_inputf(:,tt) = inhib_cum_inputf;
    spike_mat_inhibf(:, tt) = inhib_spikes_2f;
    
    % Updating the spikes for the next times step
    inhib_spikesf = inhib_spikes_2f; 
    excit_spikesf = excit_spikes_2f;
         
  
                     
end