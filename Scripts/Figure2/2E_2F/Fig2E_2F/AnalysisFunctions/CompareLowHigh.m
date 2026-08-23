%% ============================================================
% Plot the low vs high input drive states
% Produces:
%   1. pooled low/high normalization-drive histogram
%   2. empirical cumulative distributions
%   3. position-specific low/high histograms
%
% Requires from GainModulation.m:
%   normalization_drive, low_idx, high_idx, pos_bin, N_POS_BINS
% ============================================================

low_values = normalization_drive(low_idx);
high_values = normalization_drive(high_idx);

low_values = low_values(isfinite(low_values));
high_values = high_values(isfinite(high_values));

if isempty(low_values) || isempty(high_values)
    error('Low- or high-state normalization-drive data are empty.');
end

low_mean = mean(low_values, 'omitnan');
high_mean = mean(high_values, 'omitnan');
mean_difference = high_mean - low_mean;

%% Pooled Cohen's d
n_low = numel(low_values);
n_high = numel(high_values);
var_low = var(low_values, 0, 'omitnan');
var_high = var(high_values, 0, 'omitnan');

if n_low > 1 && n_high > 1
    pooled_sd = sqrt(((n_low - 1) * var_low + (n_high - 1) * var_high) / ...
                     (n_low + n_high - 2));
else
    pooled_sd = NaN;
end

if isfinite(pooled_sd) && pooled_sd > 0
    cohens_d = mean_difference / pooled_sd;
else
    cohens_d = NaN;
end

%% Histogram overlap coefficient
all_values = [low_values(:); high_values(:)];
data_min = min(all_values);
data_max = max(all_values);

if data_min == data_max
    hist_edges = [data_min - 0.5, data_max + 0.5];
else
    hist_edges = linspace(data_min, data_max, 51);
end

low_probability = histcounts(low_values, hist_edges, 'Normalization','probability');
high_probability = histcounts(high_values, hist_edges, 'Normalization','probability');
overlap_coefficient = sum(min(low_probability, high_probability));

%% Figure 1: pooled histogram
figure;
histogram(low_values, hist_edges, 'Normalization','probability', ...
    'DisplayStyle','stairs', 'LineWidth',2);
hold on
histogram(high_values, hist_edges, 'Normalization','probability', ...
    'DisplayStyle','stairs', 'LineWidth',2);
xline(low_mean, '--', 'LineWidth',1.5);
xline(high_mean, '--', 'LineWidth',1.5);
xlabel('Normalization-pool input drive (mV/s)');
ylabel('Probability');
title(sprintf('Pooled low vs high normalization drive: d = %.2f, overlap = %.2f', ...
    cohens_d, overlap_coefficient));
legend('Low-state observations','High-state observations','Low mean','High mean', ...
    'Location','best');
box off

%% Figure 2: empirical cumulative distributions
figure;
[f_low_ecdf, x_low_ecdf] = ecdf(low_values);
[f_high_ecdf, x_high_ecdf] = ecdf(high_values);
stairs(x_low_ecdf, f_low_ecdf, 'LineWidth',2);
hold on
stairs(x_high_ecdf, f_high_ecdf, 'LineWidth',2);
xlabel('Normalization-pool input drive (mV/s)');
ylabel('Cumulative probability');
title('Empirical cumulative distributions');
legend('Low state','High state','Location','best');
box off

%% Figure 3: position-specific distributions
figure;
tiledlayout(ceil(N_POS_BINS / 2), 2, 'TileSpacing','compact', 'Padding','compact');

for p = 1:N_POS_BINS
    position_mask = repmat(pos_bin == p, size(normalization_drive,1), 1);
    low_p = normalization_drive(low_idx & position_mask);
    high_p = normalization_drive(high_idx & position_mask);

    low_p = low_p(isfinite(low_p));
    high_p = high_p(isfinite(high_p));

    nexttile
    if isempty(low_p) || isempty(high_p)
        title(sprintf('Position %d: insufficient data', p));
        axis off
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

    histogram(low_p, edges_p, 'Normalization','probability', ...
        'DisplayStyle','stairs', 'LineWidth',1.5);
    hold on
    histogram(high_p, edges_p, 'Normalization','probability', ...
        'DisplayStyle','stairs', 'LineWidth',1.5);

    position_low_mean = mean(low_p, 'omitnan');
    position_high_mean = mean(high_p, 'omitnan');
    position_difference = position_high_mean - position_low_mean;

    n_low_p = numel(low_p);
    n_high_p = numel(high_p);

    if n_low_p > 1 && n_high_p > 1
        pooled_sd_p = sqrt(((n_low_p - 1) * var(low_p,0,'omitnan') + ...
                            (n_high_p - 1) * var(high_p,0,'omitnan')) / ...
                           (n_low_p + n_high_p - 2));
    else
        pooled_sd_p = NaN;
    end

    if isfinite(pooled_sd_p) && pooled_sd_p > 0
        position_cohens_d = position_difference / pooled_sd_p;
    else
        position_cohens_d = NaN;
    end

    xline(position_low_mean, '--', 'LineWidth',1);
    xline(position_high_mean, '--', 'LineWidth',1);
    xlabel('Drive (mV/s)');
    ylabel('Probability');
    title(sprintf('Position %d: \\Delta = %.2f, d = %.2f', ...
        p, position_difference, position_cohens_d));

    if p == 1
        legend('Low state','High state','Low mean','High mean','Location','best');
    end
    box off
end
