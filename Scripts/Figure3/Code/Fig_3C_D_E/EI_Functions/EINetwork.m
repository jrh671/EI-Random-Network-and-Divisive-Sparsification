if SlowTm == 1
    %% Run network simulation
        for tt = 1:length(total_time_vec)
            pos = positions(tt);
            dt = dt_vec(tt);

            if SimInputs == 1
                input_prob_firing_vec = zeros(1, n_input);
                for ii = 1:n_input
                    rfs = pf_cell{ii};
                    distance_vec = min([abs(pos - rfs); (n_pos) - abs(pos - rfs)]);
                    input_prob_firing_vec(ii) = sum(exp(-distance_vec.^2 * pf_width));
                end

                input_prob_firing_vec(input_prob_firing_vec < 0) = 0; 
                input_prob_firing_vec = (input_prob_firing_vec + alpha_P * pinknoise(n_input))* pf_rate * dt;
                input_prob_firing_vec(input_prob_firing_vec < 0) = 0; 
                coin_flips = rand(1, n_input);
                input_neurons_fired = find(coin_flips < input_prob_firing_vec & ...
                    (total_time_vec(tt) - input_most_recent_fire_times_vec) > input_refract_length);

                input_most_recent_fire_times_vec(input_neurons_fired) = total_time_vec(tt);
                input_spikes = zeros(1, n_input);
                input_spikes(input_neurons_fired) = 1;
                InputSpikes(:, tt) = input_spikes;
                clear input_spikes;
            end

            input_spikes = InputSpikes(:,tt)';

            % Updating excitatory cell membrane potential
            excit_cum_input = excit_cum_input * exp(- dt) + input_spikes * W_inputE ...
                + excit_spikes * W_EE - inhib_spikes * W_IE + alpha_Pe * pinknoise(n_excit);     
            excit_cum_input((total_time_vec(tt) - excit_most_recent_fire_times_vec) < excit_refract_length) = 0; % setting voltage equal to 0 for any neurons still in their absolute refractory period
            
            % Updating excitatory spikes 
            excit_most_recent_fire_times_vec(excit_cum_input >= I_F_thresh) = total_time_vec(tt);
            excit_spikes_2(excit_cum_input >= I_F_thresh) = 1;  % storing current timestep spikes
            excit_spikes_2(excit_cum_input < I_F_thresh) = 0;
            excit_cum_input(excit_cum_input >= I_F_thresh) = 0; % setting the input back to 0 for the neurons that spiked        
            VmE(:, tt) = excit_cum_input;
            spike_mat_excit(:, tt) = excit_spikes_2;

            % Updating inhibitory cell membrane potential
            inhib_cum_input = inhib_cum_input * exp(-dt) ... 
                + excit_spikes * W_EI + alpha_Pi * pinknoise(n_inhib);     
            inhib_cum_input((total_time_vec(tt) - inhib_most_recent_fire_times_vec) < inhib_refract_length) = 0;
            
            % Updating inhibitory spikes
            inhib_most_recent_fire_times_vec(inhib_cum_input >= I_F_thresh) = total_time_vec(tt);
            inhib_spikes_2(inhib_cum_input >= I_F_thresh) = 1;
            inhib_spikes_2(inhib_cum_input < I_F_thresh) = 0;
            inhib_cum_input(inhib_cum_input >= I_F_thresh) = 0;     
            VmI(:,tt) = inhib_cum_input;
            spike_mat_inhib(:, tt) = inhib_spikes_2;


            % Updating the spikes for the next times step
            clear excit_spikes
            clear inhib_spikes
            inhib_spikes = inhib_spikes_2; 
            excit_spikes = excit_spikes_2;
        end
else

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
end
end
