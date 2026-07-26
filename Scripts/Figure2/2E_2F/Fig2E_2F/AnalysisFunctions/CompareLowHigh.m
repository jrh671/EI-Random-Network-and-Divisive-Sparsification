%% ============================================================
% Input drive by position: low vs high normalization state
% Units: mean excitatory input drive in mV/s
%
% Note:
% input_drive_sec is summed mV drive over window_ms.
% Therefore mV/window -> mV/s is divide by window_ms.
% ============================================================

InputDrive_low  = nan(n_out,n_pos_bins);
InputDrive_high = nan(n_out,n_pos_bins);

if ~neuron_specific_split

    for p = 1:n_pos_bins

        idx_low  = low_idx  & pos_bin == p;
        idx_high = high_idx & pos_bin == p;

        if any(idx_low)
            InputDrive_low(:,p) = mean(input_drive_sec(:,idx_low),2,'omitnan');
        end

        if any(idx_high)
            InputDrive_high(:,p) = mean(input_drive_sec(:,idx_high),2,'omitnan');
        end

    end

else

    % Neuron-specific split, but measure the normalization-pool drive
    % as the input drive of all OTHER neurons, excluding neuron n.
    sum_drive = sum(input_drive_sec, 1, 'omitnan');
    count_drive = sum(isfinite(input_drive_sec), 1);

    input_drive_other_sec = nan(n_out, n_win);

    for n = 1:n_out
        other_sum = sum_drive - input_drive_sec(n,:);
        other_count = count_drive - isfinite(input_drive_sec(n,:));

        input_drive_other_sec(n,:) = other_sum ./ other_count;
        input_drive_other_sec(n, other_count <= 0) = NaN;
    end

    for n = 1:n_out

        for p = 1:n_pos_bins

            idx_low  = low_idx(n,:)  & pos_bin == p;
            idx_high = high_idx(n,:) & pos_bin == p;

            if any(idx_low)
                InputDrive_low(n,p) = mean(input_drive_other_sec(n,idx_low),'omitnan');
            end

            if any(idx_high)
                InputDrive_high(n,p) = mean(input_drive_other_sec(n,idx_high),'omitnan');
            end

        end

    end

end

mean_low_input_by_pos  = mean(InputDrive_low,1,'omitnan');
mean_high_input_by_pos = mean(InputDrive_high,1,'omitnan');

fprintf('\n');
fprintf('Position   LowInputDrive(mV/s)   HighInputDrive(mV/s)   Difference(mV/s)\n');

for p = 1:n_pos_bins
    fprintf('%3d        %20.6f   %21.6f   %18.6f\n', ...
        p, ...
        mean_low_input_by_pos(p), ...
        mean_high_input_by_pos(p), ...
        mean_high_input_by_pos(p) - mean_low_input_by_pos(p));
end


'High Input Drive (mV/s)'
mean(mean_high_input_by_pos)

'Low Input Drive (mV/s)'
mean(mean_low_input_by_pos)

'Median Difference High vs Low Input Drive (mV/s)'
median(mean_high_input_by_pos - mean_low_input_by_pos)


% DISTRIBUTION CHECK: LOW vs HIGH NORMALIZATION-POOL DRIVE
%
% Run after the main analysis and CompareLowHigh.
%
% Produces only:
%   1. Pooled low/high histogram
%   2. Kernel-density curves
%   3. Empirical cumulative distributions
%   4. Position-specific low/high histograms
% ============================================================

%% ------------------------------------------------------------
% 1. Reconstruct the normalization signal used for the split
% ------------------------------------------------------------

if ~neuron_specific_split

    % Shared population normalization pool
    normalization_drive = norm_pool;   % 1 x n_win

else

    % Leave-one-out normalization pool:
    % for neuron n, average drive of all other neurons
    sum_drive = sum(input_drive_sec, 1, 'omitnan');
    count_drive = sum(isfinite(input_drive_sec), 1);

    normalization_drive = nan(n_out, n_win);

    for n = 1:n_out

        other_sum = sum_drive - input_drive_sec(n,:);
        other_count = count_drive - isfinite(input_drive_sec(n,:));

        normalization_drive(n,:) = other_sum ./ other_count;
        normalization_drive(n, other_count <= 0) = NaN;

    end

end

%% ------------------------------------------------------------
% 2. Collect pooled low- and high-state observations
% ------------------------------------------------------------

low_values = normalization_drive(low_idx);
high_values = normalization_drive(high_idx);

low_values = low_values(isfinite(low_values));
high_values = high_values(isfinite(high_values));

if isempty(low_values) || isempty(high_values)
    error('Low- or high-state normalization-drive data are empty.');
end

fprintf('\n============================================================\n');
fprintf('LOW vs HIGH NORMALIZATION-DRIVE DISTRIBUTIONS\n');
fprintf('============================================================\n');

fprintf('Number of pooled low observations:  %d\n', numel(low_values));
fprintf('Number of pooled high observations: %d\n', numel(high_values));

low_mean = mean(low_values, 'omitnan');
low_median = median(low_values, 'omitnan');
low_sd = std(low_values, 'omitnan');

high_mean = mean(high_values, 'omitnan');
high_median = median(high_values, 'omitnan');
high_sd = std(high_values, 'omitnan');

mean_difference = high_mean - low_mean;

fprintf('\nLow drive:\n');
fprintf('  Mean   = %.4f mV/s\n', low_mean);
fprintf('  Median = %.4f mV/s\n', low_median);
fprintf('  SD     = %.4f mV/s\n', low_sd);

fprintf('\nHigh drive:\n');
fprintf('  Mean   = %.4f mV/s\n', high_mean);
fprintf('  Median = %.4f mV/s\n', high_median);
fprintf('  SD     = %.4f mV/s\n', high_sd);

fprintf('\nRaw mean difference = %.4f mV/s\n', mean_difference);

%% ------------------------------------------------------------
% 3. Calculate pooled Cohen's d
% ------------------------------------------------------------

n_low = numel(low_values);
n_high = numel(high_values);

var_low = var(low_values, 0, 'omitnan');
var_high = var(high_values, 0, 'omitnan');

if n_low > 1 && n_high > 1

    pooled_sd = sqrt( ...
        ((n_low - 1) * var_low + (n_high - 1) * var_high) / ...
        (n_low + n_high - 2));

else

    pooled_sd = NaN;

end

if isfinite(pooled_sd) && pooled_sd > 0
    cohens_d = mean_difference / pooled_sd;
else
    cohens_d = NaN;
end

fprintf('Pooled Cohen''s d = %.4f\n', cohens_d);

%% ------------------------------------------------------------
% 4. Calculate histogram overlap coefficient
%
% overlap = 1: complete overlap
% overlap = 0: no overlap
% ------------------------------------------------------------

all_values = [low_values(:); high_values(:)];

data_min = min(all_values);
data_max = max(all_values);

if data_min == data_max
    hist_edges = [data_min - 0.5, data_max + 0.5];
else
    n_bins = 50;
    hist_edges = linspace(data_min, data_max, n_bins + 1);
end

low_probability = histcounts( ...
    low_values, ...
    hist_edges, ...
    'Normalization', 'probability');

high_probability = histcounts( ...
    high_values, ...
    hist_edges, ...
    'Normalization', 'probability');

overlap_coefficient = sum(min(low_probability, high_probability));

fprintf('Histogram overlap coefficient = %.4f\n', ...
    overlap_coefficient);

%% ------------------------------------------------------------
% Figure 1. Pooled histogram
% ------------------------------------------------------------

figure;

histogram( ...
    low_values, ...
    hist_edges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

hold on;

histogram( ...
    high_values, ...
    hist_edges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

xline(low_mean, '--', 'LineWidth', 1.5);
xline(high_mean, '--', 'LineWidth', 1.5);

xlabel('Normalization-pool input drive (mV/s)');
ylabel('Probability');

title(sprintf( ...
    'Pooled low vs high normalization drive: d = %.2f, overlap = %.2f', ...
    cohens_d, overlap_coefficient));

legend( ...
    'Low-state observations', ...
    'High-state observations', ...
    'Low mean', ...
    'High mean', ...
    'Location', 'best');

box off;
hold off;



%% ------------------------------------------------------------
% Figure 2. Empirical cumulative distribution functions
% ------------------------------------------------------------

figure;

[f_low_ecdf, x_low_ecdf] = ecdf(low_values);
[f_high_ecdf, x_high_ecdf] = ecdf(high_values);

stairs(x_low_ecdf, f_low_ecdf, 'LineWidth', 2);

hold on;

stairs(x_high_ecdf, f_high_ecdf, 'LineWidth', 2);

xlabel('Normalization-pool input drive (mV/s)');
ylabel('Cumulative probability');
title('Empirical cumulative distributions');

legend( ...
    'Low state', ...
    'High state', ...
    'Location', 'best');

box off;
hold off;

%% ------------------------------------------------------------
% Figure 3. Position-specific distributions
% ------------------------------------------------------------

figure;

tiledlayout( ...
    ceil(n_pos_bins / 2), ...
    2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for p = 1:n_pos_bins

    if ~neuron_specific_split

        % Shared split: low_idx and high_idx are 1 x n_win
        idx_low_p = low_idx & (pos_bin == p);
        idx_high_p = high_idx & (pos_bin == p);

        low_p = normalization_drive(idx_low_p);
        high_p = normalization_drive(idx_high_p);

    else

        % Neuron-specific split: low_idx and high_idx are n_out x n_win
        position_mask = repmat(pos_bin == p, n_out, 1);

        idx_low_p = low_idx & position_mask;
        idx_high_p = high_idx & position_mask;

        low_p = normalization_drive(idx_low_p);
        high_p = normalization_drive(idx_high_p);

    end

    low_p = low_p(isfinite(low_p));
    high_p = high_p(isfinite(high_p));

    nexttile;

    if isempty(low_p) || isempty(high_p)

        title(sprintf('Position %d: insufficient data', p));
        axis off;
        continue

    end

    combined_p = [low_p(:); high_p(:)];

    position_min = min(combined_p);
    position_max = max(combined_p);

    if position_min == position_max
        edges_p = [position_min - 0.5, position_max + 0.5];
    else
        edges_p = linspace(position_min, position_max, 31);
    end

    histogram( ...
        low_p, ...
        edges_p, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 1.5);

    hold on;

    histogram( ...
        high_p, ...
        edges_p, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 1.5);

    position_low_mean = mean(low_p, 'omitnan');
    position_high_mean = mean(high_p, 'omitnan');
    position_difference = position_high_mean - position_low_mean;

    n_low_p = numel(low_p);
    n_high_p = numel(high_p);

    if n_low_p > 1 && n_high_p > 1

        pooled_sd_p = sqrt( ...
            ((n_low_p - 1) * var(low_p, 0, 'omitnan') + ...
             (n_high_p - 1) * var(high_p, 0, 'omitnan')) / ...
            (n_low_p + n_high_p - 2));

    else

        pooled_sd_p = NaN;

    end

    if isfinite(pooled_sd_p) && pooled_sd_p > 0
        position_cohens_d = position_difference / pooled_sd_p;
    else
        position_cohens_d = NaN;
    end

    xline(position_low_mean, '--', 'LineWidth', 1);
    xline(position_high_mean, '--', 'LineWidth', 1);

    xlabel('Drive (mV/s)');
    ylabel('Probability');

    title(sprintf( ...
        'Position %d: \\Delta = %.2f, d = %.2f', ...
        p, position_difference, position_cohens_d));

    if p == 1

        legend( ...
            'Low state', ...
            'High state', ...
            'Low mean', ...
            'High mean', ...
            'Location', 'best');

    end

    box off;
    hold off;

end