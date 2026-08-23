%% Thresholded-rate place-field overdispersion
clear SavePasses SavePDF

%% Parameters
chosen_neurons = 1:100;

activity_threshold = 2;
distance_threshold = 3;
minimum_pass_bins = 3;
minimum_expected_events = 0.25;
minimum_expected_variance = 0.10;

%% Add activity noise

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

%% Analysis settings

settings = struct( ...
    'chosen_neurons',chosen_neurons, ...
    'activity_threshold',activity_threshold, ...
    'distance_threshold',distance_threshold, ...
    'minimum_pass_bins',minimum_pass_bins, ...
    'minimum_expected_events',minimum_expected_events, ...
    'minimum_expected_variance',minimum_expected_variance);

%% Pass-wise standardized deviations

WithoutPlasticity = analyze_overdispersion_passes( ...
    RateMap1Noisy,FiringRates1Noisy,Traj1,settings);

WithPlasticity = analyze_overdispersion_passes( ...
    RateMap2Noisy,FiringRates2Noisy,Traj1,settings);

%% Plot

dataset_names = {
    'Without plasticity'
    sprintf('With plasticity; noise SD = %.3g',activity_noise_std)
};

dataset_colors = {'b','r'};

combined_z = [
    WithoutPlasticity
    WithPlasticity
];

finite_z = combined_z(isfinite(combined_z));

if isempty(finite_z)
    error('No valid passes were found.');
end

display_limit = max(3,prctile(abs(finite_z),99));
bin_edges = linspace(-display_limit,display_limit,41);

figure;
hold on;

plot_results = {WithoutPlasticity,WithPlasticity};
legend_entries = {};

for dataset_idx = 1:2

    pass_z = plot_results{dataset_idx};
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
    SavePasses{dataset_idx} = pass_z;

end

x_values = linspace(-display_limit,display_limit,500);
SavePDF = normpdf(x_values,0,1)';

plot(x_values,SavePDF,'k','LineWidth',2);

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


%% Local function: pass-wise overdispersion

function all_pass_z = analyze_overdispersion_passes(FRAve,FR,Traj,settings)

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

    for neuron_idx = neurons_to_analyze

        spatial_map = FRAve(:,:,neuron_idx);
        valid_rate_values = spatial_map(isfinite(spatial_map));

        if isempty(valid_rate_values) || max(valid_rate_values) <= 0
            continue
        end

        peak_value = max(valid_rate_values);
        [peak_rows,peak_cols] = find(spatial_map == peak_value);

        peak_coord = [mean(peak_rows) mean(peak_cols)];

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
                success_map(row_idx,col_idx) + number_of_successes;

            trial_map(row_idx,col_idx) = ...
                trial_map(row_idx,col_idx) + number_of_trials;

            time_successes(time_idx) = number_of_successes;
            time_trials(time_idx) = number_of_trials;

        end

        for pass_number = valid_pass_numbers'

            pass_indices = field_indices( ...
                pass_start_positions(pass_number): ...
                pass_end_positions(pass_number));

            observed_events = 0;
            expected_events = 0;
            expected_variance = 0;
            usable_time_bins = 0;

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

                remaining_successes = ...
                    success_map(row_idx,col_idx) - ...
                    sum(time_successes(same_bin_indices));

                remaining_trials = ...
                    trial_map(row_idx,col_idx) - ...
                    sum(time_trials(same_bin_indices));

                if remaining_trials <= 0
                    continue
                end

                event_probability = ...
                    remaining_successes / remaining_trials;

                event_probability = min(max(event_probability,eps),1-eps);

                observed_events = ...
                    observed_events + time_successes(time_idx);

                expected_events = ...
                    expected_events + ...
                    number_of_trials * event_probability;

                expected_variance = ...
                    expected_variance + ...
                    number_of_trials * event_probability * ...
                    (1-event_probability);

                usable_time_bins = usable_time_bins + 1;

            end

            if usable_time_bins == 0 || ...
               expected_events < settings.minimum_expected_events || ...
               expected_variance < settings.minimum_expected_variance
                continue
            end

            all_pass_z(end+1,1) = ...
                (observed_events - expected_events) / ...
                sqrt(expected_variance); %#ok<AGROW>

        end
    end
end


%% Local function: add Gaussian noise

function noisy_data = add_gaussian_activity_noise( ...
    input_data,noise_std,clip_at_zero)

    if isnumeric(input_data)

        noisy_data = input_data;
        finite_mask = isfinite(input_data);

        noisy_data(finite_mask) = ...
            input_data(finite_mask) + ...
            noise_std .* randn(nnz(finite_mask),1);

        if clip_at_zero
            noisy_data(finite_mask) = max(noisy_data(finite_mask),0);
        end

    elseif iscell(input_data)

        noisy_data = input_data;

        for cell_idx = 1:numel(input_data)
            noisy_data{cell_idx} = ...
                add_gaussian_activity_noise( ...
                    input_data{cell_idx},noise_std,clip_at_zero);
        end

    else

        noisy_data = input_data;

    end
end
