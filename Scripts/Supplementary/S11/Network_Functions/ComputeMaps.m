% Initialize temporal rates array
temporal_rates_CA3 = zeros(N, T);

for W = 1%[30,80]%chosen_neurons_CA3(end-1)
    chosen_neurons_CA3ave = Choose2NeuronsCA3;%[W, chosen_neurons_CA3(end)];

    % Initialize arrays to store the cumulative activity at each location for chosen neurons
    cumulative_activity_map_CA3 = zeros(P, P, length(chosen_neurons_CA3ave));
    occupancy_map = zeros(P, P);

for t = 1:T
    firing_rates_CA3 = firing_rate_snapshots_CA3{t};

    % Get the current position in the trajectory
    pos = round(trajectory(t, :)); % Round to the nearest integer to use as indices

    if all(pos > 0) && all(pos <= P) % Ensure the position is within bounds
        occupancy_map(pos(1), pos(2)) = occupancy_map(pos(1), pos(2)) + 1;

        % Update the cumulative activity map and store temporal rates
        for i = 1:N
            rate_at_pos = firing_rates_CA3{i}(pos(1), pos(2));
            
            % Update cumulative activity map if neuron is in chosen_neurons_CA3ave
            if ismember(i, chosen_neurons_CA3ave)
                neuron_index = find(chosen_neurons_CA3ave == i);
                cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) = ...
                    cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) + rate_at_pos;
            end

            % % Store the firing rate for the current time point
            temporal_rates_CA3(i, t) = rate_at_pos;
        end
    end
end

% Normalize the activity maps by the occupancy map
occupancy_map(occupancy_map == 0) = 1; % To avoid division by zero
for n = 1:length(chosen_neurons_CA3ave)
    cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) ./ occupancy_map;
end

% Normalize the activity maps by the max value of each neuron
for n = 1:length(chosen_neurons_CA3ave)
    max_val = max(cumulative_activity_map_CA3(:, :, n), [], 'all');
    if max_val > 0
        cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) / max_val;
    end
end



    % Plot the average activity maps
    % Plot for chosen neurons in CA3

  if PlotMaps==1
      figure;
    for n = 1:2
        subplot(1, 2, n);
        imagesc(cumulative_activity_map_CA3(:, :, n), [0 1]);
        colorbar;
        title(['CA3 Neuron ', num2str(chosen_neurons_CA3ave(n))]);
        set(gca, 'YDir', 'normal'); % Make y-axis upright
        xlabel('X Position');
        ylabel('Y Position');
        axis equal
        axis tight

    end

    pause(0.1)
  end
end   
pause(1)
