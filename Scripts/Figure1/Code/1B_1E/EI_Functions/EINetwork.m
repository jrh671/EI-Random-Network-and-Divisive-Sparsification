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
