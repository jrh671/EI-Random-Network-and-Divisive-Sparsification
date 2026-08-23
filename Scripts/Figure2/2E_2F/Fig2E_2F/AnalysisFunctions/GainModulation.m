%% ============================================================
% LOW-vs-HIGH GAIN ANALYSIS
%
%   - All Neurons Used
%   - remove the first 30 s (first lap)
%   - 1 s analysis windows
%   - 30 original position bins grouped into 6 groups of 5 positions
%   - leave-one-out input normalization signal
%   - within-position median split into low/high normalization states
%   - rate-defined neuron subset used for the Vm analysis
%
% Required variables in workspace:
%   InputSpikes, W_inputE, VmE, spike_mat_excit, positions,
%   neuron_indices, UnitAdjInput, dt, C
% ============================================================

%% Analysis Parameters
TOP_N = 500; %All Neurons Analyzed
WINDOW_MS = 1000;
FIRST_LAP_SEC = 30;
N_RAW_POS_BINS = 30;
N_POS_BINS = 6;
RAW_BINS_PER_GROUP = 5;
MIN_VALID_POSITIONS = 5;
MIN_VM_TUNING_VAR = 1e-4;

%% Select neurons
selected_neurons = neuron_indices(1:min(TOP_N, numel(neuron_indices)));
selected_neurons = selected_neurons(:)';

W_inputE_sel = W_inputE(:, selected_neurons);
VmE_sel = VmE(selected_neurons, :);
spike_sel = spike_mat_excit(selected_neurons, :);
positions = positions(:)';

%% Validate dimensions
T_original = size(InputSpikes, 2);

if size(VmE_sel, 2) ~= T_original
    error('VmE must have the same time dimension as InputSpikes.');
end
if size(spike_sel, 2) ~= T_original
    error('spike_mat_excit must have the same time dimension as InputSpikes.');
end
if numel(positions) ~= T_original
    error('positions must have the same length as neural activity.');
end

%% Remove first lap
bin_ms = dt * 1000;
samples_to_skip = round(FIRST_LAP_SEC * 1000 / bin_ms);

if T_original <= samples_to_skip
    error('Recording is shorter than or equal to the 30 s first lap.');
end

InputSpikes_use = InputSpikes(:, samples_to_skip+1:end);
VmE_sel_use = VmE_sel(:, samples_to_skip+1:end);
spike_sel_use = spike_sel(:, samples_to_skip+1:end);
positions_use = positions(samples_to_skip+1:end);

%% Compute weighted excitatory input drive
input_drive = (W_inputE_sel' * InputSpikes_use) * UnitAdjInput * dt / C;
n_out = size(input_drive, 1);

%% Bin all signals into 1 s windows
bins_per_window = WINDOW_MS / bin_ms;
if abs(bins_per_window - round(bins_per_window)) > eps
    error('The 1 s analysis window must contain an integer number of simulation bins.');
end
bins_per_window = round(bins_per_window);

T = size(InputSpikes_use, 2);
n_win = floor(T / bins_per_window);
T_use = n_win * bins_per_window;

if n_win < 1
    error('No complete 1 s analysis windows remain after first-lap removal.');
end

input_drive = input_drive(:, 1:T_use);
VmE_sel_use = VmE_sel_use(:, 1:T_use);
spike_sel_use = spike_sel_use(:, 1:T_use);
positions_use = positions_use(1:T_use);

VmE_centered = VmE_sel_use - mean(VmE_sel_use, 2, 'omitnan');

input_drive_rs = reshape(input_drive, n_out, bins_per_window, n_win);
VmE_centered_rs = reshape(VmE_centered, n_out, bins_per_window, n_win);
spike_rs = reshape(spike_sel_use, n_out, bins_per_window, n_win);
pos_rs = reshape(positions_use, bins_per_window, n_win);

input_drive_sec = squeeze(sum(input_drive_rs, 2));
Vm_sec_centered = squeeze(mean(VmE_centered_rs, 2, 'omitnan'));
Rate_sec = squeeze(sum(spike_rs, 2));  % spikes/s because each window is 1 s

%% Collapse 30 raw position bins into 6 groups
position_sec = mean(pos_rs, 1, 'omitnan');
pos_bin_raw = round(position_sec);
pos_bin_raw(pos_bin_raw < 1 | pos_bin_raw > N_RAW_POS_BINS) = NaN;

pos_bin = ceil(pos_bin_raw / RAW_BINS_PER_GROUP);
pos_bin(pos_bin < 1 | pos_bin > N_POS_BINS) = NaN;

%% Leave-one-out normalization signal
% For neuron n, normalization is the mean input drive of all other
% analyzed neurons in the same 1 s window.
sum_drive = sum(input_drive_sec, 1, 'omitnan');
count_drive = sum(isfinite(input_drive_sec), 1);

split_signal = nan(n_out, n_win);
for n = 1:n_out
    other_sum = sum_drive - input_drive_sec(n,:);
    other_count = count_drive - isfinite(input_drive_sec(n,:));

    split_signal(n,:) = other_sum ./ other_count;
    split_signal(n, other_count <= 0) = NaN;
end

% For when COMBINE=false
normalization_drive = split_signal;
neuron_specific_split = true;

%% Within-position median split for each neuron
low_idx = false(n_out, n_win);
high_idx = false(n_out, n_win);

for n = 1:n_out
    for p = 1:N_POS_BINS
        idx_p = (pos_bin == p) & isfinite(split_signal(n,:));
        if sum(idx_p) < 2
            continue
        end

        med_p = median(split_signal(n,idx_p), 'omitnan');
        low_idx(n,idx_p) = split_signal(n,idx_p) <= med_p;
        high_idx(n,idx_p) = split_signal(n,idx_p) > med_p;
    end
end

%% Compute low/high position tuning curves
Vm_low = nan(n_out, N_POS_BINS);
Vm_high = nan(n_out, N_POS_BINS);
Rate_low = nan(n_out, N_POS_BINS);
Rate_high = nan(n_out, N_POS_BINS);

for n = 1:n_out
    for p = 1:N_POS_BINS
        idx_low = low_idx(n,:) & (pos_bin == p);
        idx_high = high_idx(n,:) & (pos_bin == p);

        if any(idx_low)
            Vm_low(n,p) = mean(Vm_sec_centered(n,idx_low), 'omitnan');
            Rate_low(n,p) = mean(Rate_sec(n,idx_low), 'omitnan');
        end

        if any(idx_high)
            Vm_high(n,p) = mean(Vm_sec_centered(n,idx_high), 'omitnan');
            Rate_high(n,p) = mean(Rate_sec(n,idx_high), 'omitnan');
        end
    end
end

%% Define the rate-valid neuron subset and calculate Rayleigh tuning
[rate_included, rayleigh_R] = select_rate_valid_neurons( ...
    Rate_low, Rate_high, MIN_VALID_POSITIONS);

%% Apply the same neuron subset to Vm
results_Vm = analyze_vm_subset( ...
    Vm_low, Vm_high, rate_included, MIN_VALID_POSITIONS, MIN_VM_TUNING_VAR);
results_Vm.rayleigh_R = rayleigh_R;

% Kept only as a compact record of the rate selection used for Vm.
results_Rate = table(rate_included, rayleigh_R, ...
    'VariableNames', {'included','rayleigh_R'});

%% ============================================================
% Local functions
% ============================================================

function [included, rayleigh_R] = select_rate_valid_neurons(low_mat, high_mat, min_valid_positions)

    n_neurons = size(low_mat, 1);
    N_POS_BINS = size(low_mat, 2);

    included = false(n_neurons, 1);
    rayleigh_R = nan(n_neurons, 1);

    theta = linspace(0, 2*pi, N_POS_BINS + 1);
    theta = theta(1:end-1)';

    for n = 1:n_neurons
        low = low_mat(n,:)';
        high = high_mat(n,:)';

        tuning = mean([low high], 2, 'omitnan');
        good_tuning = isfinite(tuning);

        if sum(good_tuning) < min_valid_positions
            continue
        end

        tuning_use = tuning(good_tuning);
        theta_use = theta(good_tuning);

        if sum(tuning_use) <= 0
            continue
        end

        rayleigh_R(n) = abs(sum(tuning_use .* exp(1i * theta_use))) / sum(tuning_use);

        good = isfinite(low) & isfinite(high);
        low = low(good);
        high = high(good);

        if numel(low) < min_valid_positions
            continue
        end
        if var(low, 'omitnan') == 0 || var(high, 'omitnan') == 0
            continue
        end

        included(n) = true;
    end
end

function results = analyze_vm_subset(low_mat, high_mat, neuron_mask, min_valid_positions, min_tuning_var)

    n_neurons = size(low_mat, 1);

    alpha = nan(n_neurons, 1);
    beta = nan(n_neurons, 1);
    R2_affine = nan(n_neurons, 1);
    included = false(n_neurons, 1);

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

        X = [low ones(size(low))];
        coeff = X \ high;
        high_hat = X * coeff;

        alpha(n) = coeff(1);
        beta(n) = coeff(2);
        R2_affine(n) = 1 - sum((high - high_hat).^2) / ...
            sum((high - mean(high)).^2);
        included(n) = true;
    end

    results = table(alpha, beta, R2_affine, included);
end
