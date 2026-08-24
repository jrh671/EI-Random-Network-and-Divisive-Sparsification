temporal_rates_CA3 = zeros(N, T);

for W = 1
    chosen_neurons_CA3ave = Choose2NeuronsCA3;

    cumulative_activity_map_CA3 = zeros(P, P, length(chosen_neurons_CA3ave));
    occupancy_map = zeros(P, P);

for t = 1:T
    firing_rates_CA3 = firing_rate_snapshots_CA3{t};

    % Get the current position in the trajectory
    pos = round(trajectory(t, :)); 

    if all(pos > 0) && all(pos <= P) 
        occupancy_map(pos(1), pos(2)) = occupancy_map(pos(1), pos(2)) + 1;

        % Update the cumulative activity map 
        for i = 1:N
            rate_at_pos = firing_rates_CA3{i}(pos(1), pos(2));
            
            if ismember(i, chosen_neurons_CA3ave)
                neuron_index = find(chosen_neurons_CA3ave == i);
                cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) = ...
                    cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) + rate_at_pos;
            end

            temporal_rates_CA3(i, t) = rate_at_pos;
        end
    end
end

% Normalize the activity maps
occupancy_map(occupancy_map == 0) = 1; 
for n = 1:length(chosen_neurons_CA3ave)
    cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) ./ occupancy_map;
end

for n = 1:length(chosen_neurons_CA3ave)
    max_val = max(cumulative_activity_map_CA3(:, :, n), [], 'all');
    if max_val > 0
        cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) / max_val;
    end
end




if Epoch==2
  if PlotMaps==1
      figure;
    for n = 1:2
        subplot(1, 2, n);
        imagesc(cumulative_activity_map_CA3(:, :, n), [0 1]);
        SaveRateMaps{n}=cumulative_activity_map_CA3(:, :, n);
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
end

pause(1)
