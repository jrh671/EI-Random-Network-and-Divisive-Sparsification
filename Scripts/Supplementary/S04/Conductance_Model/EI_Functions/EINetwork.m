%% Run conductance-based network simulation
%
% Assumptions:
% - total_time_vec and dt_vec are in seconds.
% - Vm variables are in mV.
% - E_AMPA, E_GABA, VLeakE, VLeakI, VrE, VrI, I_F_threshE, I_F_threshI
%   are all in mV.
% - tau_AMPA, tau_GABA, and refractory lengths are all in seconds.
% - gL, W_inputE, W_EE, W_IE, W_EI, W_II, W_inputI, and C are in compatible
%   effective units.
%
% SimInputs == 1:
%   Uses normal circular place-field input.
%
% SimInputs == 0:
%   Uses controlled continuously varying input rate pf_rate_t(tt).
%   This is useful for testing divisive normalization with known input drive.
%
% Additional outputs stored by this version:
%
% InputSpikes          input-neuron spikes
% spike_mat_excit      excitatory-neuron spikes
% spike_mat_inhib      inhibitory-neuron spikes
%
% gInputE_history      input -> E conductance
% gEE_history          E -> E conductance
% gExcE                total excitatory conductance onto E cells
% gInhE                inhibitory conductance onto E cells
% gTotE                total conductance onto E cells
%
% VmE_pre              E-cell voltage before the membrane update
% VmE                  E-cell voltage after the update, before spike reset
% VinfE                instantaneous E-cell equilibrium voltage
% tau_eff_E            effective E-cell membrane time constant
%
% IInputE              input -> E current
% IEE                  E -> E current
% IExcE                total excitatory current onto E cells
% IInhE                inhibitory current onto E cells
% ILeakE               leak current onto E cells
% INetE                total current onto E cells
%
% Current sign convention:
%
%   I = g .* (E_reversal - V)
%
% Positive current drives membrane voltage upward.
% Negative current drives membrane voltage downward.

%% ------------------------------------------------------------
% Realistic synaptic reversal potentials, mV
%% ------------------------------------------------------------
E_AMPA = 0;
E_GABA = -70;

%% ------------------------------------------------------------
% Synaptic conductance decay constants, seconds
%% ------------------------------------------------------------
tau_AMPA = 0.005;
tau_GABA = 0.010;

%% ------------------------------------------------------------
% Number of simulation time points
%% ------------------------------------------------------------
n_time = length(total_time_vec);

%% ------------------------------------------------------------
% Controlled input settings for SimInputs == 0
%% ------------------------------------------------------------
T_total_sec = n_pos * n_laps;   % since each position is 1 sec
T_steps_expected = round(T_total_sec / mean(dt_vec));

% Controlled rate range in Hz
pf_rate_min = 0;
pf_rate_max = pf_rate;

% Slow continuous modulation across the whole simulation
% One full ramp up/down over the full recording.
controlled_phase = linspace(0, 2*pi, n_time);

% Continuously varying pf_rate_t in Hz
pf_rate_t = pf_rate_min + ...
    (pf_rate_max - pf_rate_min) * ...
    (0.5 + 0.5 * sin(controlled_phase));

% Save this for later analysis
InputRate_t = pf_rate_t;

%% ------------------------------------------------------------
% Conductance state variables
%% ------------------------------------------------------------
g_inputE = zeros(1, n_excit);   % input -> E
g_inputI = zeros(1, n_inhib);   % input -> I

g_EE     = zeros(1, n_excit);   % E -> E
g_IE     = zeros(1, n_excit);   % I -> E

g_EI     = zeros(1, n_inhib);   % E -> I
g_II     = zeros(1, n_inhib);   % I -> I

%% ------------------------------------------------------------
% Preallocate existing saved variables
%% ------------------------------------------------------------
InputSpikes     = zeros(n_input, n_time);
spike_mat_excit = zeros(n_excit, n_time);
spike_mat_inhib = zeros(n_inhib, n_time);

VmE = zeros(n_excit, n_time);
VmI = zeros(n_inhib, n_time);

% Existing name retained for compatibility.
% Note: I_pool stores inhibitory conductance, not inhibitory current.
I_pool = zeros(n_excit, n_time);

tau_eff_E = zeros(n_excit, n_time);

%% ------------------------------------------------------------
% Preallocate additional E-cell conductance outputs
%% ------------------------------------------------------------

% Individual excitatory conductance components
gInputE_history = zeros(n_excit, n_time);
gEE_history     = zeros(n_excit, n_time);

% Total excitatory conductance
gExcE = zeros(n_excit, n_time);

% Inhibitory conductance
gInhE = zeros(n_excit, n_time);

% Total membrane conductance
gTotE = zeros(n_excit, n_time);

%% ------------------------------------------------------------
% Preallocate additional E-cell voltage outputs
%% ------------------------------------------------------------

% Voltage at the start of the membrane update
VmE_pre = zeros(n_excit, n_time);

% Instantaneous equilibrium voltage
VinfE = zeros(n_excit, n_time);

%% ------------------------------------------------------------
% Preallocate E-cell current outputs
%
% Currents are evaluated using VmE_pre, the membrane voltage
% immediately before the conductance-based voltage update.
%% ------------------------------------------------------------

% Excitatory current components
IInputE = zeros(n_excit, n_time);
IEE     = zeros(n_excit, n_time);

% Total excitatory current
IExcE = zeros(n_excit, n_time);

% Inhibitory current
IInhE = zeros(n_excit, n_time);

% Leak current
ILeakE = zeros(n_excit, n_time);

% Net current
INetE = zeros(n_excit, n_time);

%% ------------------------------------------------------------
% Make sure previous spike vectors exist
%% ------------------------------------------------------------
if ~exist('excit_spikes', 'var') || isempty(excit_spikes)
    excit_spikes = zeros(1, n_excit);
end

if ~exist('inhib_spikes', 'var') || isempty(inhib_spikes)
    inhib_spikes = zeros(1, n_inhib);
end


pf_width = 1/pf_width_^2;
%% ------------------------------------------------------------
% Main simulation loop
%% ------------------------------------------------------------
for tt = 1:n_time

    pos = positions(tt);
    dt = dt_vec(tt);

    %% ------------------------------------------------------------
    % Input spikes
    %% ------------------------------------------------------------

    if SimInputs == 1

        %% --------------------------------------------------------
        % Normal place-field input mode
        %% --------------------------------------------------------
        input_prob_firing_vec = zeros(1, n_input);

        for ii = 1:n_input
            rfs = pf_cell{ii};

            distance_vec = min([ ...
                abs(pos - rfs); ...
                n_pos - abs(pos - rfs)]);

            input_prob_firing_vec(ii) = ...
                sum(exp(-distance_vec.^2 * pf_width));
        end

        input_prob_firing_vec(input_prob_firing_vec < 0) = 0;

        input_prob_firing_vec = ...
            (input_prob_firing_vec + ...
             alpha_P * pinknoise(n_input)) * pf_rate * dt;

        input_prob_firing_vec(input_prob_firing_vec < 0) = 0;

    else

        %% --------------------------------------------------------
        % Controlled continuous-rate input mode
        %
        % pf_rate_t(tt) is the only thing changing continuously.
        % Every input neuron gets the same global rate modulation.
        %% --------------------------------------------------------
        input_prob_firing_vec = ...
            pf_rate_t(tt) * dt * ones(1, n_input);

        % Keep valid probability range
        input_prob_firing_vec(input_prob_firing_vec < 0) = 0;
        input_prob_firing_vec(input_prob_firing_vec > 1) = 1;

    end

    %% ------------------------------------------------------------
    % Bernoulli input spike generation with refractory
    %% ------------------------------------------------------------
    coin_flips = rand(1, n_input);

    input_neurons_fired = find( ...
        coin_flips < input_prob_firing_vec & ...
        (total_time_vec(tt) - input_most_recent_fire_times_vec) ...
            > input_refract_length);

    input_most_recent_fire_times_vec(input_neurons_fired) = ...
        total_time_vec(tt);

    input_spikes = zeros(1, n_input);
    input_spikes(input_neurons_fired) = 1;

    InputSpikes(:, tt) = input_spikes(:);

    %% ------------------------------------------------------------
    % Decay synaptic conductances
    %% ------------------------------------------------------------
    g_inputE = g_inputE .* exp(-dt / tau_AMPA);
    g_inputI = g_inputI .* exp(-dt / tau_AMPA);

    g_EE     = g_EE     .* exp(-dt / tau_AMPA);
    g_EI     = g_EI     .* exp(-dt / tau_AMPA);

    g_IE     = g_IE     .* exp(-dt / tau_GABA);
    g_II     = g_II     .* exp(-dt / tau_GABA);

    %% ------------------------------------------------------------
    % Add spike-triggered conductance increments
    %% ------------------------------------------------------------

    % Feedforward input conductance onto E cells
    g_inputE = ...
        g_inputE + input_spikes * W_inputE * UnitAdjInputE;

    % Feedforward input conductance onto I cells
    g_inputI = ...
        g_inputI + input_spikes * W_inputI * UnitAdjInputI;

    % Recurrent / feedback conductances
    g_EE = ...
        g_EE + excit_spikes * W_EE;

    g_IE = ...
        g_IE + inhib_spikes * W_IE * UnitAdjInhibIE;

    % Retain original output variable
    I_pool(:, tt) = g_IE(:);

    g_EI = ...
        g_EI + excit_spikes * W_EI * UnitAdjInhibEI;

    g_II = ...
        g_II + inhib_spikes * W_II;

    %% ------------------------------------------------------------
    % Refractory masks
    %% ------------------------------------------------------------
    excit_not_refractory = ...
        (total_time_vec(tt) - excit_most_recent_fire_times_vec) ...
        >= excit_refract_length;

    inhib_not_refractory = ...
        (total_time_vec(tt) - inhib_most_recent_fire_times_vec) ...
        >= inhib_refract_length;

    %% ------------------------------------------------------------
    % Exact conductance-based update: excitatory population
    %% ------------------------------------------------------------

    % Save voltage before the membrane update.
    %
    % This is the voltage used to calculate the currents that drive
    % the membrane during this timestep.
    VmE_pre(:, tt) = excit_cum_input(:);

    % Total conductance onto each E neuron
    g_tot_E = gL + g_inputE + g_EE + g_IE;

    % Effective membrane time constant
    tau_eff_E(:, tt) = C ./ g_tot_E(:);

    % Instantaneous equilibrium voltage
    V_inf_E = ...
        (gL .* VLeakE + ...
         (g_inputE + g_EE) .* E_AMPA + ...
         g_IE .* E_GABA) ./ g_tot_E;

    %% ------------------------------------------------------------
    % Save E-cell conductance time series
    %% ------------------------------------------------------------

    % Separate excitatory conductance components
    gInputE_history(:, tt) = g_inputE(:);
    gEE_history(:, tt)     = g_EE(:);

    % Total excitatory synaptic conductance
    gExcE(:, tt) = (g_inputE + g_EE).';

    % Inhibitory synaptic conductance
    gInhE(:, tt) = g_IE(:);

    % Total conductance, including leak
    gTotE(:, tt) = g_tot_E(:);

    % Equilibrium voltage
    VinfE(:, tt) = V_inf_E(:);

    %% ------------------------------------------------------------
    % Calculate and save E-cell currents
    %
    % Sign convention:
    %
    %     I = g * (E_reversal - V)
    %
    % Currents are calculated from the pre-update membrane voltage.
    %% ------------------------------------------------------------
    V_before_E = excit_cum_input;

    % Feedforward excitatory current
    I_input_E_now = ...
        g_inputE .* (E_AMPA - V_before_E);

    % Recurrent excitatory current
    I_EE_now = ...
        g_EE .* (E_AMPA - V_before_E);

    % Total excitatory current
    I_exc_E_now = ...
        (g_inputE + g_EE) .* (E_AMPA - V_before_E);

    % Inhibitory current
    I_inh_E_now = ...
        g_IE .* (E_GABA - V_before_E);

    % Leak current
    I_leak_E_now = ...
        gL .* (VLeakE - V_before_E);

    % Net membrane current
    I_net_E_now = ...
        I_leak_E_now + I_exc_E_now + I_inh_E_now;

    % Save current time series
    IInputE(:, tt) = I_input_E_now(:);
    IEE(:, tt)     = I_EE_now(:);
    IExcE(:, tt)   = I_exc_E_now(:);
    IInhE(:, tt)   = I_inh_E_now(:);
    ILeakE(:, tt)  = I_leak_E_now(:);
    INetE(:, tt)   = I_net_E_now(:);

    %% ------------------------------------------------------------
    % Update excitatory membrane voltage
    %% ------------------------------------------------------------
    excit_cum_input(excit_not_refractory) = ...
        V_inf_E(excit_not_refractory) + ...
        (excit_cum_input(excit_not_refractory) - ...
         V_inf_E(excit_not_refractory)) .* ...
        exp( ...
            -(g_tot_E(excit_not_refractory) ./ C) * dt);

    noise_E = alpha_P * pinknoise(n_excit);
    noise_E(~excit_not_refractory) = 0;

    excit_cum_input(excit_not_refractory) = ...
        excit_cum_input(excit_not_refractory) + ...
        noise_E(excit_not_refractory);

    excit_cum_input(~excit_not_refractory) = VrE;

    %% ------------------------------------------------------------
    % Update excitatory spikes
    %% ------------------------------------------------------------
    excit_spikes_2 = zeros(1, n_excit);

    excit_spikes_2( ...
        (excit_cum_input >= I_F_threshE) & ...
        excit_not_refractory) = 1;

    excit_most_recent_fire_times_vec(excit_spikes_2 == 1) = ...
        total_time_vec(tt);

    % Save voltage after integration and before spike reset.
    VmE(:, tt) = excit_cum_input(:);

    % Reset cells that spiked
    excit_cum_input(excit_spikes_2 == 1) = VrE;

    %% ------------------------------------------------------------
    % Exact conductance-based update: inhibitory population
    %% ------------------------------------------------------------
    g_tot_I = gL + g_inputI + g_EI + g_II;

    V_inf_I = ...
        (gL .* VLeakI + ...
         g_inputI .* E_AMPA + ...
         g_EI .* E_AMPA + ...
         g_II .* E_GABA) ./ g_tot_I;

    inhib_cum_input(inhib_not_refractory) = ...
        V_inf_I(inhib_not_refractory) + ...
        (inhib_cum_input(inhib_not_refractory) - ...
         V_inf_I(inhib_not_refractory)) .* ...
        exp( ...
            -(g_tot_I(inhib_not_refractory) ./ C) * dt);

    noise_I = alpha_P * pinknoise(n_inhib);
    noise_I(~inhib_not_refractory) = 0;

    inhib_cum_input(inhib_not_refractory) = ...
        inhib_cum_input(inhib_not_refractory) + ...
        noise_I(inhib_not_refractory);

    inhib_cum_input(~inhib_not_refractory) = VrI;

    %% ------------------------------------------------------------
    % Update inhibitory spikes
    %% ------------------------------------------------------------
    inhib_spikes_2 = zeros(1, n_inhib);

    inhib_spikes_2( ...
        (inhib_cum_input >= I_F_threshI) & ...
        inhib_not_refractory) = 1;

    inhib_most_recent_fire_times_vec(inhib_spikes_2 == 1) = ...
        total_time_vec(tt);

    VmI(:, tt) = inhib_cum_input(:);

    inhib_cum_input(inhib_spikes_2 == 1) = VrI;

    %% ------------------------------------------------------------
    % Save spike variables using existing names
    %% ------------------------------------------------------------
    spike_mat_excit(:, tt) = excit_spikes_2(:);
    spike_mat_inhib(:, tt) = inhib_spikes_2(:);

    %% ------------------------------------------------------------
    % Advance spikes to next timestep
    %% ------------------------------------------------------------
    excit_spikes = excit_spikes_2;
    inhib_spikes = inhib_spikes_2;

    %% ------------------------------------------------------------
    % Progress display
    %% ------------------------------------------------------------
    steps_per_second = round(1 / dt);

    if steps_per_second > 0 && ...
            mod(tt, steps_per_second) == 0 && ...
            tt > steps_per_second

        win = tt - steps_per_second + 1 : tt;

        total_spikes = sum(spike_mat_excit(:, win), 'all');

        n_active_cells = ...
            sum(any(spike_mat_excit(:, win), 2));

        % Avoid division by zero when no E cells fired.
        if n_active_cells > 0
            spikes_per_active_cell = ...
                total_spikes / n_active_cells;
        else
            spikes_per_active_cell = 0;
        end

        fprintf( ...
            't = %.1f s | spikes/active cell = %.2f | active cells = %d\n', ...
            total_time_vec(tt), ...
            spikes_per_active_cell, ...
            n_active_cells);
    end
end