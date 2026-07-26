%% ============================================================
% LOW-vs-HIGH NORMALIZATION POOL ANALYSIS
% VmE + firing rates
%
% Required variables in workspace:
%   InputSpikes
%   W_inputE
%   VmE
%   spike_mat_excit
%   positions
%   neuron_indices
%   UnitAdjInput
%   dt
%   C
%
% Modes:
%   NormPoolMode = 'population':
%       Original analysis. Low/high state is defined by the population-average
%       input drive across selected neurons. This tests a shared normalization pool.
%
%   NormPoolMode = 'neuron_specific':
%       Control analysis. Low/high state is defined separately for each neuron
%       using that neuron's own input_drive_sec(n,:). This tests explicit
%       neuron-specific drive modulation, not a shared population pool.
%
%   PositionSpecific = 1:
%       Original position-specific tuning-curve analysis.
%       Low/high normalization split is done within each position bin.
%       Then low/high position tuning curves are compared using gain/offset fits.
%
%   PositionSpecific = 0:
%       True global low/high analysis.
%       Low/high normalization split is done globally.
%       Position is NOT used to compute tuning curves.
%       Each neuron gets one global low mean and one global high mean.
% ============================================================


PositionSpecific = 1;

% 'population' or 'neuron_specific'
% population: original shared normalization pool
% neuron_specific: each neuron gets its own low/high split based on its own drive
NormPoolMode = 'neuron_specific';

%% -------------------------
% Parameters
% -------------------------

bin_ms = dt*1000;
window_ms = 1000;
bins_per_window = window_ms / bin_ms;
window_sec = window_ms / 1000;

n_pos_bins = 30;
pos_group_size = 5;              % 3 adjacent original bins per grouped bin
n_pos_bins_original = n_pos_bins;
n_pos_bins = n_pos_bins_original / pos_group_size;

top_n = 500;

% Optional first-lap removal
remove_first_lap = true;
first_lap_sec = n_pos_bins_original;

% Low/high normalization split
low_prctile = 50;
high_prctile = 50;
split_mode = 'median';   % 'median' or 'quantile'

% Position-specific tuning parameters
min_samples_per_pos = 1;
min_valid_positions = 5;

% Rate activity filters for PositionSpecific == 1
min_mean_rate = 0;
min_peak_rate = 0;
min_rate_range = 0;
min_rate_var = 0;

% Tuning filter for PositionSpecific == 1
rayleigh_thresh = 0;

% Vm reliability filter for PositionSpecific == 1
min_vm_tuning_var = 1e-4;

% Global-analysis filters for PositionSpecific == 0
min_global_mean_rate = 0;
min_global_peak_rate = 0;
min_global_vm_abs_diff = 0;

%% -------------------------
% Select top neurons
% -------------------------

selected_neurons = neuron_indices(1:min(top_n, numel(neuron_indices)));
selected_neurons = selected_neurons(:)';

W_inputE_sel = W_inputE(:, selected_neurons);
VmE_sel = VmE(selected_neurons, :);
spike_sel = spike_mat_excit(selected_neurons, :);

%% -------------------------
% Format/check
% -------------------------

positions = positions(:)';

T_original = size(InputSpikes, 2);

if size(VmE_sel, 2) ~= T_original
    error('VmE must have same time dimension as InputSpikes.');
end

if size(spike_sel, 2) ~= T_original
    error('spike_mat_excit must have same time dimension as InputSpikes.');
end

if numel(positions) ~= T_original
    error('positions must have same length as neural activity.');
end

%% -------------------------
% Create analysis copies
% -------------------------

InputSpikes_use = InputSpikes;
VmE_sel_use = VmE_sel;
spike_sel_use = spike_sel;
positions_use = positions;

if remove_first_lap

    samples_to_skip = round(first_lap_sec * 1000 / bin_ms);

    if T_original <= samples_to_skip
        error('Recording is shorter than or equal to first_lap_sec.');
    end

    InputSpikes_use = InputSpikes_use(:, samples_to_skip+1:end);
    VmE_sel_use = VmE_sel_use(:, samples_to_skip+1:end);
    spike_sel_use = spike_sel_use(:, samples_to_skip+1:end);
    positions_use = positions_use(samples_to_skip+1:end);

    fprintf('Skipped first %.1f s = %d samples.\n', ...
        first_lap_sec, samples_to_skip);

else

    fprintf('First lap was NOT skipped.\n');

end

T = size(InputSpikes_use, 2);

%% ============================================================
% 1. Compute weighted input drive
% ============================================================

% input_drive = W_inputE_sel' * InputSpikes_use;
input_drive = (W_inputE_sel' * InputSpikes_use) * UnitAdjInput * dt / C;
n_out = size(input_drive, 1);

%% ============================================================
% 2. Bin into analysis windows
% ============================================================

n_win = floor(T / bins_per_window);
T_use = n_win * bins_per_window;

input_drive = input_drive(:, 1:T_use);
VmE_sel_use = VmE_sel_use(:, 1:T_use);
spike_sel_use = spike_sel_use(:, 1:T_use);
positions_use = positions_use(1:T_use);

VmE_centered = VmE_sel_use - mean(VmE_sel_use, 2, 'omitnan');

input_drive_rs = reshape(input_drive, n_out, bins_per_window, n_win);
VmE_centered_rs = reshape(VmE_centered, n_out, bins_per_window, n_win);
VmE_raw_rs = reshape(VmE_sel_use, n_out, bins_per_window, n_win);
spike_rs = reshape(spike_sel_use, n_out, bins_per_window, n_win);
pos_rs = reshape(positions_use, bins_per_window, n_win);

input_drive_sec = squeeze(sum(input_drive_rs, 2));

% Centered Vm is used only for position-specific tuning curves, preserving
% the original PositionSpecific == 1 analysis behavior.
Vm_sec_centered = squeeze(mean(VmE_centered_rs, 2, 'omitnan'));

% Raw Vm is used for global low/high means, because centered Vm forces
% low and high means to balance around zero when the split is even.
Vm_sec_raw = squeeze(mean(VmE_raw_rs, 2, 'omitnan'));

Rate_sec = squeeze(sum(spike_rs, 2)) / window_sec;

% position_sec = mean(pos_rs, 1, 'omitnan');
% pos_bin = round(position_sec);
% pos_bin(pos_bin < 1 | pos_bin > n_pos_bins) = NaN;

position_sec = mean(pos_rs, 1, 'omitnan');
pos_bin_raw = round(position_sec);

pos_bin_raw(pos_bin_raw < 1 | pos_bin_raw > n_pos_bins_original) = NaN;

% Group adjacent position bins:
% original bins 1-3   -> grouped bin 1
% original bins 4-6   -> grouped bin 2
% ...
% original bins 28-30 -> grouped bin 10
pos_bin = ceil(pos_bin_raw / pos_group_size);
pos_bin(pos_bin < 1 | pos_bin > n_pos_bins) = NaN;


%% ============================================================
% 3. Define normalization pool and low/high windows
% ============================================================

% Population-level normalization pool. This is always kept as norm_pool so
% downstream regression/plotting scripts that expect a 1 x n_win vector keep
% working exactly as before.
norm_pool = mean(input_drive_sec, 1, 'omitnan');

% The signal used only for splitting low/high windows can either be the
% population pool or each neuron's own input drive.
switch NormPoolMode

    case 'population'

        split_signal = norm_pool;              % 1 x n_win
        low_idx = false(1, n_win);             % shared windows
        high_idx = false(1, n_win);
        neuron_specific_split = false;

    case 'neuron_specific'

    % Leave-one-out normalization pool:
    % for neuron n, split by the mean input drive of all OTHER neurons.
    sum_drive = sum(input_drive_sec, 1, 'omitnan');
    count_drive = sum(isfinite(input_drive_sec), 1);

    split_signal = nan(n_out, n_win);

    for n = 1:n_out
        other_sum = sum_drive - input_drive_sec(n,:);
        other_count = count_drive - isfinite(input_drive_sec(n,:));

        split_signal(n,:) = other_sum ./ other_count;
        split_signal(n, other_count <= 0) = NaN;
    end

    low_idx = false(n_out, n_win);
    high_idx = false(n_out, n_win);
    neuron_specific_split = true;

    otherwise

        error('NormPoolMode must be either population or neuron_specific.');

end

if PositionSpecific == 1

    %% ------------------------------------------------------------
    % Position-stratified low/high split
    % ------------------------------------------------------------

    if ~neuron_specific_split

        % Original behavior: one population normalization pool defines the
        % same low/high windows for all neurons, within each position bin.
        for p = 1:n_pos_bins

            idx_p = isfinite(pos_bin) & pos_bin == p & isfinite(split_signal);

            if sum(idx_p) < 2
                continue
            end

            Np = split_signal(idx_p);

            switch split_mode

                case 'median'

                    med_p = median(Np, 'omitnan');

                    low_idx(idx_p)  = split_signal(idx_p) <= med_p;
                    high_idx(idx_p) = split_signal(idx_p) > med_p;

                case 'quantile'

                    low_thresh_p = prctile(Np, low_prctile);
                    high_thresh_p = prctile(Np, high_prctile);

                    low_idx(idx_p)  = split_signal(idx_p) <= low_thresh_p;
                    high_idx(idx_p) = split_signal(idx_p) >= high_thresh_p;

                otherwise

                    error('split_mode must be either median or quantile.');

            end

        end

        low_idx = low_idx & isfinite(pos_bin);
        high_idx = high_idx & isfinite(pos_bin);

        fprintf('Analysis mode: POSITION-SPECIFIC tuning curves.\n');
        fprintf('Split mode: POPULATION normalization pool.\n');
        fprintf('Low-N windows: %d\n', sum(low_idx));
        fprintf('High-N windows: %d\n', sum(high_idx));

    else

        % Neuron-specific behavior: each neuron gets its own low/high windows
        % based on its own input drive, within each position bin.
        for n = 1:n_out

            for p = 1:n_pos_bins

                idx_p = isfinite(pos_bin) & pos_bin == p & isfinite(split_signal(n,:));

                if sum(idx_p) < 2
                    continue
                end

                Np = split_signal(n,idx_p);

                switch split_mode

                    case 'median'

                        med_p = median(Np, 'omitnan');

                        low_idx(n,idx_p)  = split_signal(n,idx_p) <= med_p;
                        high_idx(n,idx_p) = split_signal(n,idx_p) > med_p;

                    case 'quantile'

                        low_thresh_p = prctile(Np, low_prctile);
                        high_thresh_p = prctile(Np, high_prctile);

                        low_idx(n,idx_p)  = split_signal(n,idx_p) <= low_thresh_p;
                        high_idx(n,idx_p) = split_signal(n,idx_p) >= high_thresh_p;

                    otherwise

                        error('split_mode must be either median or quantile.');

                end

            end

            low_idx(n,:) = low_idx(n,:) & isfinite(pos_bin);
            high_idx(n,:) = high_idx(n,:) & isfinite(pos_bin);

        end

        low_counts = sum(low_idx, 2);
        high_counts = sum(high_idx, 2);

        fprintf('Analysis mode: POSITION-SPECIFIC tuning curves.\n');
        fprintf('Split mode: NEURON-SPECIFIC own input drive.\n');
        fprintf('Low-N windows per neuron: mean %.1f, median %.1f\n', ...
            mean(low_counts, 'omitnan'), median(low_counts, 'omitnan'));
        fprintf('High-N windows per neuron: mean %.1f, median %.1f\n', ...
            mean(high_counts, 'omitnan'), median(high_counts, 'omitnan'));

    end

else

    %% ------------------------------------------------------------
    % Global low/high split
    % Position is ignored for the response analysis.
    % ------------------------------------------------------------

    if ~neuron_specific_split

        % Original behavior: one population split used for all neurons.
        valid_idx = isfinite(split_signal);

        switch split_mode

            case 'median'

                med_N = median(split_signal(valid_idx), 'omitnan');

                low_idx(valid_idx)  = split_signal(valid_idx) <= med_N;
                high_idx(valid_idx) = split_signal(valid_idx) >  med_N;

            case 'quantile'

                low_thresh_N = prctile(split_signal(valid_idx), low_prctile);
                high_thresh_N = prctile(split_signal(valid_idx), high_prctile);

                low_idx(valid_idx)  = split_signal(valid_idx) <= low_thresh_N;
                high_idx(valid_idx) = split_signal(valid_idx) >= high_thresh_N;

            otherwise

                error('split_mode must be either median or quantile.');

        end

        fprintf('Analysis mode: GLOBAL low/high means; position ignored.\n');
        fprintf('Split mode: POPULATION normalization pool.\n');
        fprintf('Low-N windows: %d\n', sum(low_idx));
        fprintf('High-N windows: %d\n', sum(high_idx));

    else

        % Neuron-specific global split: each neuron is split by its own
        % input drive, ignoring position.
        for n = 1:n_out

            valid_idx = isfinite(split_signal(n,:));

            switch split_mode

                case 'median'

                    med_N = median(split_signal(n,valid_idx), 'omitnan');

                    low_idx(n,valid_idx)  = split_signal(n,valid_idx) <= med_N;
                    high_idx(n,valid_idx) = split_signal(n,valid_idx) >  med_N;

                case 'quantile'

                    low_thresh_N = prctile(split_signal(n,valid_idx), low_prctile);
                    high_thresh_N = prctile(split_signal(n,valid_idx), high_prctile);

                    low_idx(n,valid_idx)  = split_signal(n,valid_idx) <= low_thresh_N;
                    high_idx(n,valid_idx) = split_signal(n,valid_idx) >= high_thresh_N;

                otherwise

                    error('split_mode must be either median or quantile.');

            end

        end

        low_counts = sum(low_idx, 2);
        high_counts = sum(high_idx, 2);

        fprintf('Analysis mode: GLOBAL low/high means; position ignored.\n');
        fprintf('Split mode: NEURON-SPECIFIC own input drive.\n');
        fprintf('Low-N windows per neuron: mean %.1f, median %.1f\n', ...
            mean(low_counts, 'omitnan'), median(low_counts, 'omitnan'));
        fprintf('High-N windows per neuron: mean %.1f, median %.1f\n', ...
            mean(high_counts, 'omitnan'), median(high_counts, 'omitnan'));

    end

end

%% ============================================================
% 4. Run selected analysis mode
% ============================================================

if PositionSpecific == 1

    %% ============================================================
    % 4A. Compute low/high position tuning curves
    % ============================================================

    Vm_low = nan(n_out, n_pos_bins);
    Vm_high = nan(n_out, n_pos_bins);

    Rate_low = nan(n_out, n_pos_bins);
    Rate_high = nan(n_out, n_pos_bins);

    n_low_by_pos = zeros(1, n_pos_bins);
    n_high_by_pos = zeros(1, n_pos_bins);

    if ~neuron_specific_split

        % Original behavior: same low/high windows for every neuron.
        for p = 1:n_pos_bins

            idx_low = low_idx & pos_bin == p;
            idx_high = high_idx & pos_bin == p;

            n_low_by_pos(p) = sum(idx_low);
            n_high_by_pos(p) = sum(idx_high);

            if n_low_by_pos(p) >= min_samples_per_pos
                Vm_low(:,p) = mean(Vm_sec_centered(:,idx_low), 2, 'omitnan');
                Rate_low(:,p) = mean(Rate_sec(:,idx_low), 2, 'omitnan');
            end

            if n_high_by_pos(p) >= min_samples_per_pos
                Vm_high(:,p) = mean(Vm_sec_centered(:,idx_high), 2, 'omitnan');
                Rate_high(:,p) = mean(Rate_sec(:,idx_high), 2, 'omitnan');
            end

        end

    else

        % Neuron-specific behavior: each neuron uses its own low/high windows.
        n_low_by_pos_neuron = zeros(n_out, n_pos_bins);
        n_high_by_pos_neuron = zeros(n_out, n_pos_bins);

        for n = 1:n_out

            for p = 1:n_pos_bins

                idx_low = low_idx(n,:) & pos_bin == p;
                idx_high = high_idx(n,:) & pos_bin == p;

                n_low_by_pos_neuron(n,p) = sum(idx_low);
                n_high_by_pos_neuron(n,p) = sum(idx_high);

                if n_low_by_pos_neuron(n,p) >= min_samples_per_pos
                    Vm_low(n,p) = mean(Vm_sec_centered(n,idx_low), 2, 'omitnan');
                    Rate_low(n,p) = mean(Rate_sec(n,idx_low), 2, 'omitnan');
                end

                if n_high_by_pos_neuron(n,p) >= min_samples_per_pos
                    Vm_high(n,p) = mean(Vm_sec_centered(n,idx_high), 2, 'omitnan');
                    Rate_high(n,p) = mean(Rate_sec(n,idx_high), 2, 'omitnan');
                end

            end

        end

        n_low_by_pos = mean(n_low_by_pos_neuron, 1, 'omitnan');
        n_high_by_pos = mean(n_high_by_pos_neuron, 1, 'omitnan');

    end

    %% ============================================================
    % 5A. Analyze rates first: active + Rayleigh tuned neurons only
    % ============================================================

    results_Rate = analyze_gain_offset_active_rates_rayleigh( ...
        Rate_low, Rate_high, ...
        min_valid_positions, ...
        min_mean_rate, ...
        min_peak_rate, ...
        min_rate_range, ...
        min_rate_var, ...
        rayleigh_thresh);

    active_tuned_mask = results_Rate.included;

    fprintf('\nActive + Rayleigh-tuned rate neurons retained: %d / %d\n', ...
        sum(active_tuned_mask), n_out);

    %% ============================================================
    % 6A. Analyze Vm using the SAME final neuron subset
    % ============================================================

    results_Vm = analyze_gain_offset_subset( ...
        Vm_low, Vm_high, ...
        active_tuned_mask, ...
        min_valid_positions, ...
        min_vm_tuning_var);

    results_Vm.rayleigh_R = results_Rate.rayleigh_R;

    %% ============================================================
    % 7A. Summary
    % ============================================================

    Vm_summary = summarize_gain_offset(results_Vm, 'Vm_active_Rayleigh');
    Rate_summary = summarize_gain_offset(results_Rate, 'Rate_active_Rayleigh');

else

    %% ============================================================
    % 4B. Compute global low/high means
    % ============================================================

    if ~neuron_specific_split

        % Original behavior: same global low/high windows for every neuron.
        Vm_low_global = mean(Vm_sec_raw(:,low_idx), 2, 'omitnan');
        Vm_high_global = mean(Vm_sec_raw(:,high_idx), 2, 'omitnan');

        Rate_low_global = mean(Rate_sec(:,low_idx), 2, 'omitnan');
        Rate_high_global = mean(Rate_sec(:,high_idx), 2, 'omitnan');

    else

        % Neuron-specific behavior: each neuron uses its own global low/high windows.
        Vm_low_global = nan(n_out, 1);
        Vm_high_global = nan(n_out, 1);

        Rate_low_global = nan(n_out, 1);
        Rate_high_global = nan(n_out, 1);

        for n = 1:n_out

            Vm_low_global(n) = mean(Vm_sec_raw(n,low_idx(n,:)), 2, 'omitnan');
            Vm_high_global(n) = mean(Vm_sec_raw(n,high_idx(n,:)), 2, 'omitnan');

            Rate_low_global(n) = mean(Rate_sec(n,low_idx(n,:)), 2, 'omitnan');
            Rate_high_global(n) = mean(Rate_sec(n,high_idx(n,:)), 2, 'omitnan');

        end

    end

    Vm_diff_global = Vm_high_global - Vm_low_global;
    Rate_diff_global = Rate_high_global - Rate_low_global;

    Vm_ratio_global = Vm_high_global ./ Vm_low_global;
    Rate_ratio_global = Rate_high_global ./ Rate_low_global;

    Vm_ratio_global(abs(Vm_low_global) <= eps) = NaN;
    Rate_ratio_global(abs(Rate_low_global) <= eps) = NaN;

    %% ============================================================
    % 5B. Analyze global low/high differences
    % ============================================================

    results_Vm = analyze_global_low_high( ...
        Vm_low_global, ...
        Vm_high_global, ...
        min_global_vm_abs_diff, ...
        'Vm_global');

    results_Rate = analyze_global_low_high_rate( ...
        Rate_low_global, ...
        Rate_high_global, ...
        min_global_mean_rate, ...
        min_global_peak_rate, ...
        'Rate_global');

    fprintf('\nGlobal Vm neurons retained: %d / %d\n', ...
        sum(results_Vm.included), n_out);

    fprintf('Global rate neurons retained: %d / %d\n', ...
        sum(results_Rate.included), n_out);

    %% ============================================================
    % 6B. Summary
    % ============================================================

    Vm_summary = summarize_global_low_high(results_Vm, 'Vm_global');
    Rate_summary = summarize_global_low_high(results_Rate, 'Rate_global');

end

%% ============================================================
% 8. Display summaries
% ============================================================

disp(Vm_summary)
disp(Rate_summary)

CompareLowHigh;

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function results = analyze_gain_offset_active_rates_rayleigh( ...
    low_mat, high_mat, min_valid_positions, min_mean_rate, ...
    min_peak_rate, min_rate_range, min_rate_var, rayleigh_thresh)

    n_neurons = size(low_mat,1);
    n_pos_bins = size(low_mat,2);

    alpha = nan(n_neurons,1);
    beta = nan(n_neurons,1);
    R2_affine = nan(n_neurons,1);

    peak_low = nan(n_neurons,1);
    peak_high = nan(n_neurons,1);
    peak_ratio = nan(n_neurons,1);

    mean_low = nan(n_neurons,1);
    mean_high = nan(n_neurons,1);
    mean_diff = nan(n_neurons,1);

    shape_corr = nan(n_neurons,1);
    norm_shape_dist = nan(n_neurons,1);
    rayleigh_R = nan(n_neurons,1);

    included = false(n_neurons,1);

    theta = linspace(0, 2*pi, n_pos_bins+1);
    theta = theta(1:end-1)';

    for n = 1:n_neurons

        low_full = low_mat(n,:)';
        high_full = high_mat(n,:)';

        tuning_all = mean([low_full high_full], 2, 'omitnan');
        good_tuning = isfinite(tuning_all);

        if sum(good_tuning) < min_valid_positions
            continue
        end

        tuning_use = tuning_all(good_tuning);
        theta_use = theta(good_tuning);

        if sum(tuning_use) <= 0
            continue
        end

        rayleigh_R(n) = abs(sum(tuning_use .* exp(1i * theta_use))) / ...
                        sum(tuning_use);

        if rayleigh_R(n) < rayleigh_thresh
            continue
        end

        low = low_full;
        high = high_full;

        good = isfinite(low) & isfinite(high);
        low = low(good);
        high = high(good);

        if numel(low) < min_valid_positions
            continue
        end

        both = [low; high];

        mean_rate = mean(both, 'omitnan');
        peak_rate = max(both);
        rate_range = max(both) - min(both);
        rate_var = var(both, 'omitnan');

        if mean_rate < min_mean_rate
            continue
        end

        if peak_rate < min_peak_rate
            continue
        end

        if rate_range < min_rate_range
            continue
        end

        if rate_var < min_rate_var
            continue
        end

        if var(low, 'omitnan') == 0 || var(high, 'omitnan') == 0
            continue
        end

        [alpha(n), beta(n), R2_affine(n)] = affine_fit(low, high);

        peak_low(n) = max(low);
        peak_high(n) = max(high);

        if peak_low(n) > eps
            peak_ratio(n) = peak_high(n) / peak_low(n);
        end

        mean_low(n) = mean(low, 'omitnan');
        mean_high(n) = mean(high, 'omitnan');
        mean_diff(n) = mean_high(n) - mean_low(n);

        [shape_corr(n), norm_shape_dist(n)] = shape_metrics(low, high);

        included(n) = true;

    end

    results = make_results_table(alpha, beta, R2_affine, peak_low, peak_high, ...
        peak_ratio, mean_low, mean_high, mean_diff, shape_corr, ...
        norm_shape_dist, included);

    results.rayleigh_R = rayleigh_R;

end

function results = analyze_gain_offset_subset( ...
    low_mat, high_mat, neuron_mask, min_valid_positions, min_tuning_var)

    n_neurons = size(low_mat,1);

    alpha = nan(n_neurons,1);
    beta = nan(n_neurons,1);
    R2_affine = nan(n_neurons,1);

    peak_low = nan(n_neurons,1);
    peak_high = nan(n_neurons,1);
    peak_ratio = nan(n_neurons,1);

    mean_low = nan(n_neurons,1);
    mean_high = nan(n_neurons,1);
    mean_diff = nan(n_neurons,1);

    shape_corr = nan(n_neurons,1);
    norm_shape_dist = nan(n_neurons,1);

    included = false(n_neurons,1);

    for n = 1:n_neurons

        if ~neuron_mask(n)
            continue
        end

        low = low_mat(n,:)';
        high = high_mat(n,:)';

        good = isfinite(low) & isfinite(high);
        low = low(good);
        high = high(good);

        if numel(low) < min_valid_positions
            continue
        end

        if var(low, 'omitnan') < min_tuning_var
            continue
        end

        if var(high, 'omitnan') == 0
            continue
        end

        [alpha(n), beta(n), R2_affine(n)] = affine_fit(low, high);

        peak_low(n) = max(low);
        peak_high(n) = max(high);

        if abs(peak_low(n)) > eps
            peak_ratio(n) = peak_high(n) / peak_low(n);
        end

        mean_low(n) = mean(low, 'omitnan');
        mean_high(n) = mean(high, 'omitnan');
        mean_diff(n) = mean_high(n) - mean_low(n);

        [shape_corr(n), norm_shape_dist(n)] = shape_metrics(low, high);

        included(n) = true;

    end

    results = make_results_table(alpha, beta, R2_affine, peak_low, peak_high, ...
        peak_ratio, mean_low, mean_high, mean_diff, shape_corr, ...
        norm_shape_dist, included);

end

function [alpha, beta, R2] = affine_fit(low, high)

    X = [low ones(size(low))];
    coeff = X \ high;

    alpha = coeff(1);
    beta = coeff(2);

    high_hat = X * coeff;

    R2 = 1 - sum((high - high_hat).^2) / ...
        sum((high - mean(high)).^2);

end

function [shape_corr, norm_shape_dist] = shape_metrics(low, high)

    low_z = low - mean(low, 'omitnan');
    high_z = high - mean(high, 'omitnan');

    low_norm = norm(low_z);
    high_norm = norm(high_z);

    if low_norm > 0 && high_norm > 0
        low_shape = low_z / low_norm;
        high_shape = high_z / high_norm;

        shape_corr = corr(low_shape, high_shape, 'rows', 'complete');
        norm_shape_dist = norm(high_shape - low_shape);
    else
        shape_corr = NaN;
        norm_shape_dist = NaN;
    end

end

function results = make_results_table(alpha, beta, R2_affine, peak_low, ...
    peak_high, peak_ratio, mean_low, mean_high, mean_diff, ...
    shape_corr, norm_shape_dist, included)

    results = table;
    results.alpha = alpha;
    results.beta = beta;
    results.R2_affine = R2_affine;
    results.peak_low = peak_low;
    results.peak_high = peak_high;
    results.peak_ratio = peak_ratio;
    results.mean_low = mean_low;
    results.mean_high = mean_high;
    results.mean_diff = mean_diff;
    results.shape_corr = shape_corr;
    results.norm_shape_dist = norm_shape_dist;
    results.included = included;

end

function summary = summarize_gain_offset(results, label)

    valid = results.included & isfinite(results.alpha);

    summary = table;

    summary.Signal = {label};
    summary.N = sum(valid);

    summary.Mean_alpha = mean(results.alpha(valid), 'omitnan');
    summary.Median_alpha = median(results.alpha(valid), 'omitnan');

    summary.Mean_beta = mean(results.beta(valid), 'omitnan');
    summary.Median_beta = median(results.beta(valid), 'omitnan');

    summary.Mean_R2_affine = mean(results.R2_affine(valid), 'omitnan');

    summary.Mean_peak_ratio = mean(results.peak_ratio(valid), 'omitnan');
    summary.Median_peak_ratio = median(results.peak_ratio(valid), 'omitnan');

    summary.Mean_mean_diff = mean(results.mean_diff(valid), 'omitnan');
    summary.Median_mean_diff = median(results.mean_diff(valid), 'omitnan');

    summary.Mean_shape_corr = mean(results.shape_corr(valid), 'omitnan');
    summary.Median_shape_corr = median(results.shape_corr(valid), 'omitnan');

    if ismember('rayleigh_R', results.Properties.VariableNames)
        summary.Mean_Rayleigh_R = mean(results.rayleigh_R(valid), 'omitnan');
        summary.Median_Rayleigh_R = median(results.rayleigh_R(valid), 'omitnan');
    end

end

function results = analyze_global_low_high(low_vec, high_vec, min_abs_diff, label)

    low_vec = low_vec(:);
    high_vec = high_vec(:);

    good = isfinite(low_vec) & isfinite(high_vec);

    diff_vec = high_vec - low_vec;
    ratio_vec = high_vec ./ low_vec;
    ratio_vec(abs(low_vec) <= eps) = NaN;

    included = good & abs(diff_vec) >= min_abs_diff;

    results = table;
    results.Signal = repmat({label}, numel(low_vec), 1);
    results.low_mean = low_vec;
    results.high_mean = high_vec;
    results.diff = diff_vec;
    results.ratio = ratio_vec;
    results.included = included;

end

function results = analyze_global_low_high_rate(low_vec, high_vec, min_mean_rate, min_peak_rate, label)

    low_vec = low_vec(:);
    high_vec = high_vec(:);

    good = isfinite(low_vec) & isfinite(high_vec);

    diff_vec = high_vec - low_vec;
    ratio_vec = high_vec ./ low_vec;
    ratio_vec(abs(low_vec) <= eps) = NaN;

    mean_rate = mean([low_vec high_vec], 2, 'omitnan');
    peak_rate = max([low_vec high_vec], [], 2);

    included = good & mean_rate >= min_mean_rate & peak_rate >= min_peak_rate;

    results = table;
    results.Signal = repmat({label}, numel(low_vec), 1);
    results.low_mean = low_vec;
    results.high_mean = high_vec;
    results.diff = diff_vec;
    results.ratio = ratio_vec;
    results.mean_rate = mean_rate;
    results.peak_rate = peak_rate;
    results.included = included;

end

function summary = summarize_global_low_high(results, label)

    valid = results.included;

    summary = table;

    summary.Signal = {label};
    summary.N = sum(valid);

    summary.Mean_low = mean(results.low_mean(valid), 'omitnan');
    summary.Median_low = median(results.low_mean(valid), 'omitnan');

    summary.Mean_high = mean(results.high_mean(valid), 'omitnan');
    summary.Median_high = median(results.high_mean(valid), 'omitnan');

    summary.Mean_diff = mean(results.diff(valid), 'omitnan');
    summary.Median_diff = median(results.diff(valid), 'omitnan');

    summary.Mean_ratio = mean(results.ratio(valid), 'omitnan');
    summary.Median_ratio = median(results.ratio(valid), 'omitnan');

end
