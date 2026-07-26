        
if Across == 1
        %Iteration for loop
        W_EI = (Iteration-1)*0.015 * ones(n_excit, n_inhib);
        W_IE = (Iteration-1)*0.015 * ones(n_inhib, n_excit);
else
        W_EI = 0.05 * ones(n_excit, n_inhib);
        W_IE = 0.05 * ones(n_inhib, n_excit);
end

    if SlowTm == 1

        %Integrate and Fire Parameters
        I_F_thresh = 1; 
        input_refract_length = 0;
        excit_refract_length = 0.02;
        inhib_refract_length = 0.02;

        %Connection Scaling Factors
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

    else        
        %Integrate and Fire Parameters        
        I_F_threshE = -55;                           % integrate and fire threshold 
        I_F_threshI = -55;                           % integrate and fire threshold 
        
        VL = -70;
        VrE = -65;
        VrI = -65;
        VLeakE = VL*ones(1,n_excit); %Leak Potential
        VLeakI = VL*ones(1,n_inhib); %Leak Potential
        
        C=.0000005; %Capacitance u(mF) where mF = (Amp/mV) | Note S.I. F = (Amp/V)
        gL=.000025; %Conductance u(mS) where mS = (C/mV)   | Note S.I. S = (C/V)
        Tm=(C/gL); %Resulting units are in seconds
        input_refract_length = 0;                   % input neuron refractory length 
        excit_refract_length = 0.002;                % excit refractory length
        inhib_refract_length = 0.002;                % inhibitory refractory length
        
        %Connection Scaling Factors
        prob_pf_center = 0.20; 
        prob_E_recurrent_connectivity = 0.25;       % probability of any given excitatory neuron will connect to another excitatory neuron
        prob_I_E_connectivity = 0.30;               % probability of any inhibitory neuron connecting to an excitatory neuron
        prob_E_I_connectivity = 0.175;
        prob_I_I_connectivity = .50;              % probability of any I neuron connecting to an inhibitory neuron   
        
        %Connection Limits
        InputWeight=1; 
        UnitAdjInput = .00015; 
        
        IIWeight=1;
        UnitAdjInhibEI = .075; 
        UnitAdjInhibIE = .075/20;
        
        MaxWeight=0.5;
        WeightLim=MaxWeight;
        initial_weight_max_input = 0.05;            % initial weight maximum input to E
        initial_weight_max_EE = 0.05 / prob_E_recurrent_connectivity; % initial weight maximum for E to E
        initial_weight_max_EI = 0.05 / prob_I_E_connectivity;  % initial weight maximum for E to I
        initial_weight_max_IE = 0.05 / prob_E_I_connectivity; % initial weight maximum for I to E
        initial_weight_max_II = 0.05/ prob_I_I_connectivity;
        W_upper_limit = MaxWeight;                        % maximum weight we allow the E to I weights to take
        
        pf_rate = 25;                                % input neuron rate (in Hz) at center of pf
        pf_width = 1/6;                              % controls the Gaussian pf width (note it's actually the inverse) 
    end

        if Across == 2 %Input Noise
            alpha_P = (Iteration-1)*0.02;
        else
            alpha_P = 0;
        end
        
        if Across == 3 %Output Noise
            alpha_Pe = (Iteration-1)*0.02;
        else
            alpha_Pe = 0; 
        end
        
        alpha_Pi = 0; %Inhibition Noise

        %% Initialize weights and place fields
        InputStrength = 1;

        if Iteration==1
        W_inputE = InputStrength * rand(n_input, n_excit) * initial_weight_max_input;
        end
       
        W_EE = 0;
        W_II = 0;
        n_laps=n_laps+1;
