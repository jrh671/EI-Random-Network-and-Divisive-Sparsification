temporal_rates_CA3 = zeros(N, T);

Plot=0;

for Epoch=1:2

firing_rate_snapshots_CA3=Rates{Epoch};

if Epoch == 1
trajectory=Trajectory11;
elseif Epoch == 2
trajectory=Trajectory1;    
end

chosen_neurons_CA3ave = 1:1:100;   % should contain up to 100 neuron indices

if length(chosen_neurons_CA3ave) > 100
    chosen_neurons_CA3ave = chosen_neurons_CA3ave(1:100);
end

cumulative_activity_map_CA3 = zeros(P, P, length(chosen_neurons_CA3ave));
occupancy_map = zeros(P, P);

for t = 1:T
    firing_rates_CA3 = firing_rate_snapshots_CA3{t};

    % Get the current position
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

% Normalize the activity maps by the occupancy map
occupancy_map(occupancy_map == 0) = 1; % Avoid division by zero
for n = 1:length(chosen_neurons_CA3ave)
    cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) ./ occupancy_map;
end

for n = 1:length(chosen_neurons_CA3ave)
    max_val = max(cumulative_activity_map_CA3(:, :, n), [], 'all');
    if max_val > 0
        cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) / max_val;
    end
end

% Plot the average activity maps
if Epoch == 2
    if Plot == 1

        nPerFigure = 10;   
        nFiguresMax = 10;  
        totalNeurons = length(chosen_neurons_CA3ave);

        neuronsToPlot = min(totalNeurons, nPerFigure * nFiguresMax);

        SaveRateMaps = cell(1, neuronsToPlot);

        for fig_idx = 1:ceil(neuronsToPlot / nPerFigure)

            figure;

            start_idx = (fig_idx - 1) * nPerFigure + 1;
            end_idx   = min(fig_idx * nPerFigure, neuronsToPlot);

            for plot_idx = start_idx:end_idx
                local_idx = plot_idx - start_idx + 1;

                subplot(2, 5, local_idx);
                imagesc(cumulative_activity_map_CA3(:, :, plot_idx), [0 1]);
                SaveRateMaps{plot_idx} = cumulative_activity_map_CA3(:, :, plot_idx);
                colorbar;
                title(['CA3 Neuron ', num2str(chosen_neurons_CA3ave(plot_idx))]);
                set(gca, 'YDir', 'normal');
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

GoalAccumulation;

end