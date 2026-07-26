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

clear SheetOverdisp OverdispersionResults ...
      PlasticityOverdispersionNoiseCheck

%% Parameters

chosen_neurons = 1:100;

activity_threshold = 2;
distance_threshold = 3;
minimum_pass_bins = 3;
minimum_expected_events = 0;
minimum_expected_variance = 0;
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

% %% Package analysis settings
% 
% analysis_settings = struct( ...
%     'chosen_neurons',chosen_neurons, ...
%     'activity_threshold',activity_threshold, ...
%     'distance_threshold',distance_threshold, ...
%     'minimum_pass_bins',minimum_pass_bins, ...
%     'minimum_expected_events',minimum_expected_events, ...
%     'minimum_expected_variance',minimum_expected_variance, ...
%     'minimum_passes_per_neuron',minimum_passes_per_neuron, ...
%     'smoothing_alpha',smoothing_alpha, ...
%     'smoothing_beta',smoothing_beta);
% 
% %% Run the analyses
% 
% WithoutPlasticityWithNoise = analyze_overdispersion_dataset( ...
%     RateMap1Noisy,FiringRates1Noisy,Traj1, ...
%     sprintf('Without plasticity: noise SD = %.3g',activity_noise_std), ...
%     analysis_settings);
% 
% PlasticityWithNoise = analyze_overdispersion_dataset( ...
%     RateMap2Noisy,FiringRates2Noisy,Traj2, ...
%     sprintf('With plasticity: noise SD = %.3g',activity_noise_std), ...
%     analysis_settings);
% 
% % This extra no-noise run is used ONLY for the plasticity variance check.
% PlasticityWithoutNoise = analyze_overdispersion_dataset( ...
%     RateMap2,FiringRates2,Traj2, ...
%     'With plasticity: no added noise',analysis_settings);

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

end

x_values = linspace(-display_limit,display_limit,500);

plot(x_values,normpdf(x_values,0,1), ...
    'k','LineWidth',2);

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
