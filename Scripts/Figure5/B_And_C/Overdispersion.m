%% Thresholded-rate place-field overdispersion analysis
%
% Required inputs:
%   RateMap1, RateMap2         row x column x neuron
%   FiringRates1, FiringRates2 FR{time_idx}{neuron_idx}
%   Traj1, Traj2               time x 2, [row, column]
%
% This version:
%   1. Adds Gaussian noise to BOTH experimental conditions.
%   2. Uses both noisy conditions for the main overdispersion plots and
%      statistical comparison.
%   3. Separately analyzes the original, no-noise plasticity condition.
%   4. Compares the final plasticity neuron-level overdispersion variance
%      without noise versus with noise, using neurons valid in both runs.
%
% No files are written.

clear SheetOverdisp OverdispersionResults ...
      PlasticityOverdispersionNoiseCheck

%% Parameters

chosen_neurons = 1:100;

activity_threshold = 2;
distance_threshold = 3;
minimum_pass_bins = 3;
minimum_expected_events = 0.25;
minimum_expected_variance = 0.10;
minimum_passes_per_neuron = 5;

smoothing_alpha = 0;
smoothing_beta  = 0;

%% Noise settings: applied to BOTH datasets

activity_noise_std = 1;
clip_noisy_activity_at_zero = true;
activity_noise_seed = 1;

rng(activity_noise_seed,'twister');

RateMap1Noisy = add_gaussian_activity_noise( ...
    RateMap1,activity_noise_std,clip_noisy_activity_at_zero);

RateMap2Noisy = add_gaussian_activity_noise( ...
    RateMap2,activity_noise_std,clip_noisy_activity_at_zero);

FiringRates1Noisy = add_gaussian_activity_noise( ...
    FiringRates1,activity_noise_std,clip_noisy_activity_at_zero);

FiringRates2Noisy = add_gaussian_activity_noise( ...
    FiringRates2,activity_noise_std,clip_noisy_activity_at_zero);

%% Package analysis settings

analysis_settings = struct( ...
    'chosen_neurons',chosen_neurons, ...
    'activity_threshold',activity_threshold, ...
    'distance_threshold',distance_threshold, ...
    'minimum_pass_bins',minimum_pass_bins, ...
    'minimum_expected_events',minimum_expected_events, ...
    'minimum_expected_variance',minimum_expected_variance, ...
    'minimum_passes_per_neuron',minimum_passes_per_neuron, ...
    'smoothing_alpha',smoothing_alpha, ...
    'smoothing_beta',smoothing_beta);

%% Run the analyses

WithoutPlasticityWithNoise = analyze_overdispersion_dataset( ...
    RateMap1Noisy,FiringRates1Noisy,Traj1, ...
    sprintf('Without plasticity: noise SD = %.3g',activity_noise_std), ...
    analysis_settings);

PlasticityWithNoise = analyze_overdispersion_dataset( ...
    RateMap2Noisy,FiringRates2Noisy,Traj2, ...
    sprintf('With plasticity: noise SD = %.3g',activity_noise_std), ...
    analysis_settings);

% This extra no-noise run is used ONLY for the plasticity variance check.
PlasticityWithoutNoise = analyze_overdispersion_dataset( ...
    RateMap2,FiringRates2,Traj2, ...
    'With plasticity: no added noise',analysis_settings);

%% Preserve the original-style outputs
%
% The main condition comparison uses:
%   dataset 1: without plasticity, no added noise
%   dataset 2: with plasticity, with added noise

OverdispersionResults(1) = WithoutPlasticityWithNoise;
OverdispersionResults(2) = PlasticityWithNoise;

SheetOverdisp = {
    WithoutPlasticityWithNoise.AllPassZ
    PlasticityWithNoise.AllPassZ
};

%% Plot pass-wise standardized deviations

dataset_names = {
    'Without plasticity'
    sprintf('With plasticity; noise SD = %.3g',activity_noise_std)
};

dataset_colors = {'b','r'};

figure;
hold on;

combined_z = [
    WithoutPlasticityWithNoise.AllPassZ
    PlasticityWithNoise.AllPassZ
];

finite_z = combined_z(isfinite(combined_z));

if isempty(finite_z)
    error('No valid passes were found.');
end

display_limit = max(3,prctile(abs(finite_z),99));
bin_edges = linspace(-display_limit,display_limit,41);
legend_entries = {};

plot_results = {WithoutPlasticityWithNoise,PlasticityWithNoise};

for dataset_idx = 1:2

    pass_z = plot_results{dataset_idx}.AllPassZ;
    pass_z = pass_z(isfinite(pass_z));

    if isempty(pass_z)
        continue
    end

    histogram(pass_z, ...
        'Normalization','pdf', ...
        'BinEdges',bin_edges, ...
        'FaceAlpha',0.40, ...
        'FaceColor',dataset_colors{dataset_idx});

    legend_entries{end+1} = dataset_names{dataset_idx}; %#ok<SAGROW>

    SavePasses{dataset_idx}=pass_z;
end

x_values = linspace(-display_limit,display_limit,500);

plot(x_values,normpdf(x_values,0,1), ...
    'k','LineWidth',2);

SavePDF=normpdf(x_values,0,1)';

legend_entries{end+1} = 'Independent-event expectation';

xline(0,'k:','LineWidth',1);

xlabel('Pass-wise standardized event deviation');
ylabel('Probability density');

title({ ...
    'Place-field overdispersion of thresholded rate events'
    sprintf('Activity threshold = %.3g',activity_threshold)
});

legend(legend_entries,'Location','best');
box on;
hold off;

%% Plot neuron-level overdispersion for the main condition comparison

without_plasticity = ...
    WithoutPlasticityWithNoise.Neurons.Overdispersion;

with_plasticity = ...
    PlasticityWithNoise.Neurons.Overdispersion;

without_plasticity = ...
    without_plasticity(isfinite(without_plasticity));

with_plasticity = ...
    with_plasticity(isfinite(with_plasticity));

figure;
hold on;

group_values = [
    ones(size(without_plasticity))
    2*ones(size(with_plasticity))
];

overdispersion_values = [
    without_plasticity
    with_plasticity
];

if ~isempty(overdispersion_values)

    boxchart(group_values,overdispersion_values);

    yline(1,'k--', ...
        'Independent-event expectation', ...
        'LineWidth',1.5);

end

xlim([0.5 2.5]);
xticks([1 2]);
xticklabels(dataset_names);

ylabel('Overdispersion: variance of pass-wise deviations');
title('Thresholded-rate place-field overdispersion by neuron');

box on;
hold off;

%% Compare the two main conditions statistically

if ~isempty(without_plasticity) && ~isempty(with_plasticity)

    p_overdispersion = ranksum( ...
        without_plasticity,with_plasticity);

    fprintf('\nNeuron-level overdispersion comparison\n');

    fprintf(['Without plasticity + noise: median = %.3f, ' ...
             'n = %d neurons\n'], ...
        median(without_plasticity), ...
        numel(without_plasticity));

    fprintf(['With plasticity + noise: median = %.3f, ' ...
             'n = %d neurons\n'], ...
        median(with_plasticity), ...
        numel(with_plasticity));

    fprintf('Rank-sum p-value = %.6g\n',p_overdispersion);

else

    p_overdispersion = NaN;

    warning(['Insufficient valid neurons for the ' ...
             'between-condition comparison.']);

end

%% Final overdispersion variance check for plasticity only
%
% Compare the final plotted neuron-level overdispersion values from:
%   - RateMap2 / FiringRates2
%   - RateMap2Noisy / FiringRates2Noisy
%
% Neurons are matched by neuron number so both variances use the same set.

table_without_noise = PlasticityWithoutNoise.Neurons;
table_with_noise = PlasticityWithNoise.Neurons;

[common_neurons,index_without,index_with] = intersect( ...
    table_without_noise.Neuron, ...
    table_with_noise.Neuron, ...
    'stable');

matched_without_noise = ...
    table_without_noise.Overdispersion(index_without);

matched_with_noise = ...
    table_with_noise.Overdispersion(index_with);

valid_pair = ...
    isfinite(matched_without_noise) & ...
    isfinite(matched_with_noise);

common_neurons = common_neurons(valid_pair);
matched_without_noise = matched_without_noise(valid_pair);
matched_with_noise = matched_with_noise(valid_pair);

if numel(common_neurons) >= 2

    variance_without_noise = var(matched_without_noise,0);
    variance_with_noise = var(matched_with_noise,0);

    variance_growth = ...
        variance_with_noise - variance_without_noise;

    percent_variance_growth = ...
        100 * variance_growth / variance_without_noise;

    variance_ratio = ...
        variance_with_noise / variance_without_noise;

else

    variance_without_noise = NaN;
    variance_with_noise = NaN;
    variance_growth = NaN;
    percent_variance_growth = NaN;
    variance_ratio = NaN;

    warning(['Fewer than two matched plasticity neurons were valid ' ...
             'in both the no-noise and noisy analyses.']);

end

PlasticityOverdispersionNoiseCheck = struct( ...
    'NoiseStandardDeviation',activity_noise_std, ...
    'NoiseSeed',activity_noise_seed, ...
    'MatchedNeurons',common_neurons, ...
    'OverdispersionWithoutNoise',matched_without_noise, ...
    'OverdispersionWithNoise',matched_with_noise, ...
    'VarianceWithoutNoise',variance_without_noise, ...
    'VarianceWithNoise',variance_with_noise, ...
    'VarianceGrowth',variance_growth, ...
    'PercentVarianceGrowth',percent_variance_growth, ...
    'VarianceRatio',variance_ratio);

fprintf('\n');
fprintf('====================================================\n');
fprintf('FINAL PLASTICITY OVERDISPERSION VARIANCE CHECK\n');
fprintf('====================================================\n');
fprintf('Noise SD:                       %.6f\n', ...
    activity_noise_std);
fprintf('Matched neurons:                %d\n', ...
    numel(common_neurons));
fprintf('Variance without added noise:   %.6f\n', ...
    variance_without_noise);
fprintf('Variance with added noise:      %.6f\n', ...
    variance_with_noise);
fprintf('Variance growth:                %.6f\n', ...
    variance_growth);
fprintf('Percent variance growth:        %.2f%%\n', ...
    percent_variance_growth);
fprintf('Variance ratio, noisy/original: %.4f\n', ...
    variance_ratio);
fprintf('====================================================\n');



%% Local function: analyze one dataset

function result = analyze_overdispersion_dataset( ...
    FRAve,FR,Traj,dataset_name,settings)

    number_of_time_bins = min(size(Traj,1),numel(FR));
    Traj = Traj(1:number_of_time_bins,:);

    map_rows = size(FRAve,1);
    map_cols = size(FRAve,2);
    number_of_neurons = size(FRAve,3);

    neurons_to_analyze = settings.chosen_neurons( ...
        settings.chosen_neurons >= 1 & ...
        settings.chosen_neurons <= number_of_neurons);

    trajectory_rows = round(Traj(:,1));
    trajectory_cols = round(Traj(:,2));

    valid_position = ...
        isfinite(trajectory_rows) & ...
        isfinite(trajectory_cols) & ...
        trajectory_rows >= 1 & ...
        trajectory_rows <= map_rows & ...
        trajectory_cols >= 1 & ...
        trajectory_cols <= map_cols;

    all_pass_z = [];
    all_pass_neuron_ids = [];

    neuron_pass_count = zeros(number_of_neurons,1);
    neuron_mean_z = nan(number_of_neurons,1);
    neuron_overdispersion = nan(number_of_neurons,1);

    pass_results = [];

    fprintf('\nAnalyzing %s\n',dataset_name);

    for neuron_idx = neurons_to_analyze

        spatial_map = FRAve(:,:,neuron_idx);
        valid_rate_values = spatial_map(isfinite(spatial_map));

        if isempty(valid_rate_values) || max(valid_rate_values) <= 0
            continue
        end

        peak_value = max(valid_rate_values);
        [peak_rows,peak_cols] = find(spatial_map == peak_value);

        peak_coord = [
            mean(peak_rows)
            mean(peak_cols)
        ]';

        distances = nan(number_of_time_bins,1);

        distances(valid_position) = vecnorm( ...
            Traj(valid_position,:) - peak_coord,2,2);

        inside_field = ...
            valid_position & ...
            distances <= settings.distance_threshold;

        field_indices = find(inside_field);

        if isempty(field_indices)
            continue
        end

        gap_locations = find(diff(field_indices) > 1);

        pass_start_positions = [1; gap_locations + 1];
        pass_end_positions = [gap_locations; numel(field_indices)];

        pass_lengths = ...
            pass_end_positions - pass_start_positions + 1;

        valid_pass_numbers = find( ...
            pass_lengths >= settings.minimum_pass_bins);

        if isempty(valid_pass_numbers)
            continue
        end

        success_map = zeros(map_rows,map_cols);
        trial_map = zeros(map_rows,map_cols);

        time_successes = zeros(number_of_time_bins,1);
        time_trials = zeros(number_of_time_bins,1);

        for time_idx = 1:number_of_time_bins

            if ~valid_position(time_idx)
                continue
            end

            if isempty(FR{time_idx}) || ...
               neuron_idx > numel(FR{time_idx}) || ...
               isempty(FR{time_idx}{neuron_idx})
                continue
            end

            rate_values = FR{time_idx}{neuron_idx}(:);
            rate_values = rate_values(isfinite(rate_values));

            if isempty(rate_values)
                continue
            end

            number_of_trials = numel(rate_values);

            number_of_successes = sum( ...
                rate_values >= settings.activity_threshold);

            row_idx = trajectory_rows(time_idx);
            col_idx = trajectory_cols(time_idx);

            success_map(row_idx,col_idx) = ...
                success_map(row_idx,col_idx) + ...
                number_of_successes;

            trial_map(row_idx,col_idx) = ...
                trial_map(row_idx,col_idx) + ...
                number_of_trials;

            time_successes(time_idx) = number_of_successes;
            time_trials(time_idx) = number_of_trials;

        end

        neuron_pass_z = [];

        for pass_number = valid_pass_numbers'

            pass_indices = field_indices( ...
                pass_start_positions(pass_number): ...
                pass_end_positions(pass_number));

            observed_events = 0;
            expected_events = 0;
            expected_variance = 0;
            usable_time_bins = 0;
            total_trials_in_pass = 0;

            for time_idx = pass_indices'

                number_of_trials = time_trials(time_idx);

                if number_of_trials <= 0
                    continue
                end

                row_idx = trajectory_rows(time_idx);
                col_idx = trajectory_cols(time_idx);

                same_bin_in_pass = ...
                    trajectory_rows(pass_indices) == row_idx & ...
                    trajectory_cols(pass_indices) == col_idx;

                same_bin_indices = pass_indices(same_bin_in_pass);

                removed_successes = ...
                    sum(time_successes(same_bin_indices));

                removed_trials = ...
                    sum(time_trials(same_bin_indices));

                remaining_successes = ...
                    success_map(row_idx,col_idx) - removed_successes;

                remaining_trials = ...
                    trial_map(row_idx,col_idx) - removed_trials;

                probability_denominator = ...
                    remaining_trials + ...
                    settings.smoothing_alpha + ...
                    settings.smoothing_beta;

                if probability_denominator <= 0
                    continue
                end

                event_probability = ...
                    (remaining_successes + ...
                     settings.smoothing_alpha) / ...
                    probability_denominator;

                event_probability = min(max( ...
                    event_probability,eps),1-eps);

                observed_events = ...
                    observed_events + time_successes(time_idx);

                expected_events = ...
                    expected_events + ...
                    number_of_trials * event_probability;

                expected_variance = ...
                    expected_variance + ...
                    number_of_trials * ...
                    event_probability * ...
                    (1-event_probability);

                total_trials_in_pass = ...
                    total_trials_in_pass + number_of_trials;

                usable_time_bins = usable_time_bins + 1;

            end

            if usable_time_bins == 0 || ...
               expected_events < settings.minimum_expected_events || ...
               expected_variance < settings.minimum_expected_variance
                continue
            end

            pass_z = ...
                (observed_events - expected_events) / ...
                sqrt(expected_variance);

            neuron_pass_z(end+1,1) = pass_z; %#ok<AGROW>
            all_pass_z(end+1,1) = pass_z; %#ok<AGROW>
            all_pass_neuron_ids(end+1,1) = neuron_idx; %#ok<AGROW>

            pass_results(end+1,:) = [ ... %#ok<AGROW>
                neuron_idx
                pass_number
                pass_indices(1)
                pass_indices(end)
                usable_time_bins
                total_trials_in_pass
                observed_events
                expected_events
                expected_variance
                pass_z
            ]';

        end

        neuron_pass_count(neuron_idx) = numel(neuron_pass_z);

        if numel(neuron_pass_z) >= ...
                settings.minimum_passes_per_neuron

            neuron_mean_z(neuron_idx) = mean(neuron_pass_z);

            neuron_overdispersion(neuron_idx) = ...
                var(neuron_pass_z,0);

        end

    end

    pass_variable_names = {
        'Neuron'
        'PassNumber'
        'StartTimeBin'
        'EndTimeBin'
        'UsableTimeBins'
        'NumberOfRateSamples'
        'ObservedEvents'
        'ExpectedEvents'
        'ExpectedVariance'
        'StandardizedDeviation'
    };

    if isempty(pass_results)

        pass_table = array2table( ...
            zeros(0,numel(pass_variable_names)), ...
            'VariableNames',pass_variable_names);

    else

        pass_table = array2table( ...
            pass_results, ...
            'VariableNames',pass_variable_names);

    end

    valid_neuron_indices = find( ...
        neuron_pass_count >= settings.minimum_passes_per_neuron);

    neuron_table = table( ...
        valid_neuron_indices, ...
        neuron_pass_count(valid_neuron_indices), ...
        neuron_mean_z(valid_neuron_indices), ...
        neuron_overdispersion(valid_neuron_indices), ...
        'VariableNames',{
            'Neuron'
            'NumberOfPasses'
            'MeanStandardizedDeviation'
            'Overdispersion'
        });

    result = struct;
    result.Dataset = dataset_name;
    result.Passes = pass_table;
    result.Neurons = neuron_table;
    result.AllPassZ = all_pass_z;
    result.PassNeuronIDs = all_pass_neuron_ids;

    fprintf('  Valid standardized passes: %d\n', ...
        numel(all_pass_z));

    fprintf('  Neurons with at least %d passes: %d\n', ...
        settings.minimum_passes_per_neuron, ...
        height(neuron_table));

    if ~isempty(neuron_table)

        fprintf('  Median neuron overdispersion: %.3f\n', ...
            median(neuron_table.Overdispersion,'omitnan'));

    end

end

%% Local function: recursively add Gaussian noise

function noisy_data = add_gaussian_activity_noise( ...
    input_data,noise_std,clip_at_zero)

    if isnumeric(input_data)

        noisy_data = input_data;
        finite_mask = isfinite(input_data);

        noisy_data(finite_mask) = ...
            input_data(finite_mask) + ...
            noise_std .* randn(nnz(finite_mask),1);

        if clip_at_zero

            noisy_data(finite_mask) = max( ...
                noisy_data(finite_mask),0);

        end

    elseif iscell(input_data)

        noisy_data = input_data;

        for cell_idx = 1:numel(input_data)

            noisy_data{cell_idx} = ...
                add_gaussian_activity_noise( ...
                input_data{cell_idx}, ...
                noise_std,clip_at_zero);

        end

    else

        noisy_data = input_data;

    end

end
