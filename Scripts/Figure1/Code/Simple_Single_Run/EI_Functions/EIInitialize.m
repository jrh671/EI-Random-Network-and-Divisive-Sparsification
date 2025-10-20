        %Integrate and Fire Parameters
        I_F_thresh = 1; 
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

        pf_rate = 12; % Input neuron rate (in Hz) at center of pf
        pf_width = 10; % Controls Gaussian pf width pf_width = (1/Sigma^2)
        alpha_P = 0; %Noise Input
        alpha_Pe = 0; %Noise E
        alpha_Pi = 0; %Noise I

        %% Initialize weights and place fields
        InputStrength = 1;
        W_inputE = InputStrength * rand(n_input, n_excit) * initial_weight_max_input;
        W_EE = 0;
        W_EI = 0.05 * ones(n_excit, n_inhib);
        W_IE = 0.05 * ones(n_inhib, n_excit);
        W_II = 0;
