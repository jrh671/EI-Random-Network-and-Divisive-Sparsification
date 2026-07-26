
% Find the top 4 neurons with the highest average activity at the special location for CA1
chosen_neurons_CA3 = sorted_indices_CA3(1:8);
chosen_neurons_CA3(length(chosen_neurons_CA3)-1:length(chosen_neurons_CA3)) = Choose2NeuronsCA3;

ComputeMaps;

if VisualizePlasticity == 1

figure;
plot(trajectory(:, 1), trajectory(:, 2), 'o-');
if special_radius > 0
hold on;
plot(special_location(1), special_location(2), 'rx', 'MarkerSize', 12, 'LineWidth', 2);  % Mark the special location
viscircles(special_location, special_radius, 'Color', 'r'); % Draw the special location radius
hold off;
end
xlabel('X Position');
ylabel('Y Position');
title('Full Combined Trajectory');
grid on;
axis equal
axis tight
pause(1)

if Plasticity==1
% Visualization of the chosen neurons' activities evolving through time
figure('units','normalized','outerposition',[0 0 1 1])

for t = TLearning:EndTime
    for n = 1:2
        subplot(1, 2, n);
        p_t = trajectory(t, :);
        chosen_neuron = chosen_neurons_CA3(end-2+n);
        activity = firing_rate_snapshots_CA3{t}{chosen_neuron};
        imagesc(activity);
        set(gca, 'YDir', 'normal');  % Make y-axis upright
        hold on;
        plot(p_t(1), p_t(2), 'rx', 'MarkerSize', 25, 'LineWidth', 2);  % Mark the current position
        if special_radius > 0
            plot(special_location(1), special_location(2), 'go', 'MarkerSize', 8, 'LineWidth', 2);  % Mark the special location
            viscircles(special_location, special_radius, 'Color', 'g'); % Draw the special location radius
        end
        title(sprintf('Neuron %d Activity at Time %d', chosen_neuron, t));
        hold off;
        colorbar;
        caxis([0, CLimits(n)]);
        axis equal
        axis tight
    end
    pause(0.1);  % Pause to create an animation effect
end
end
end

