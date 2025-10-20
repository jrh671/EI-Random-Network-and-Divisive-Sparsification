figure;

%%
hold on;

colors = {'b', 'r'}; % Colors for Dataset 1 and Dataset 2
legend_entries = {};


for i = 1:2
    if i == 1
        FRAve = RateMap1;
        FR = FiringRates1;
        Traj = Traj1;
        dataset_name = 'Dataset 1';
    else
        FRAve = RateMap2;
        FR = FiringRates2;
        Traj = Traj2;
        dataset_name = 'Dataset 2';
    end

    % Define parameters
    chosen_neurons = 1:34;%nonZeroIndices1'; 
    distance_threshold = 3; 
    activity_threshold = 8;

    all_firing_rates = [];
    bin_edges = linspace(-1.5, 1.5, 30); % Fixed bin edges for both histograms

    for neuron_idx = chosen_neurons
        spatial_map = FRAve(:,:,neuron_idx);

        % Find all indices with max value
        max_value = max(spatial_map(:));
        peak_indices = find(spatial_map == max_value);
        
        % Randomly select one max position (or choose first occurrence)
        chosen_peak_idx = peak_indices(randi(length(peak_indices))); % Choose at random if multiple PFs
        
        % Convert linear index to subscripts
        [peak_row, peak_col] = ind2sub(size(spatial_map), chosen_peak_idx);
        peak_coord = [peak_row, peak_col];
    
        % Compute distances
        distances = vecnorm(Traj - peak_coord, 2, 2);
        pass_indices = find(distances <= distance_threshold);
    
        % Keep only indices where at least Consec consecutive time bins are within the threshold
        Consec = 0; % Set the required minimum consecutive time bins
        diff_pass = diff(pass_indices);
        split_idx = find(diff_pass > 1); % Identify gaps in continuity
        start_idx = [1; split_idx + 1];
        end_idx = [split_idx; length(pass_indices)];
    
    
        % Filter sequences that are at least N bins long
        valid_sequences = (end_idx - start_idx + 1) >= Consec;
        valid_passes = [];
        for j = find(valid_sequences)'
            valid_passes = [valid_passes; pass_indices(start_idx(j):end_idx(j))];
        end
    
            neuron_firing_rates = [];
            for idx = valid_passes'
    
                
                % Define noise scaling factor 
                noise_std = 0.2; % Uniform Noise
                
                % Add Gaussian noise
                snapshot = FR{idx}{neuron_idx}(:) + noise_std * randn(size(FR{idx}{neuron_idx}(:)));
                snapshot = snapshot(snapshot >= activity_threshold);
                neuron_firing_rates = [neuron_firing_rates; snapshot];
            end
    
            all_firing_rates = [all_firing_rates; neuron_firing_rates];
    end

    if i==2
        all_firing_rates=all_firing_rates(randperm(length(all_firing_rates),length(zscored_firing_rates)));
    end

    length(all_firing_rates)
    if ~isempty(all_firing_rates)
        zscored_firing_rates = zscore(all_firing_rates);

        histogram(zscored_firing_rates,'Normalization','pdf', ...
                  'FaceAlpha', 0.5, 'FaceColor', colors{i}, 'BinEdges', bin_edges);

       SheetOverdisp{i} = zscored_firing_rates;
        
        legend_entries{end+1} = dataset_name;
    else
        disp(['No valid firing rates above threshold for ' dataset_name '.']);
    end
end

% Overlay theoretical Gaussian curve
x_vals = linspace(-3.5, 3.5, 100);
y_vals = normpdf(x_vals, 0.5, .5);
plot(x_vals, y_vals, 'k', 'LineWidth', 2);
legend_entries{end+1} = 'Theoretical Gaussian';



% Labels and legend
xlabel('Z-scored Firing Rate');
ylabel('Density');
title('Z-scored Firing Rates from Spatial Passes (Each Neuron’s Own Peak)');
legend(legend_entries);
hold off;

xlim([-1.5,1.5]);
ylim([0,1]);
