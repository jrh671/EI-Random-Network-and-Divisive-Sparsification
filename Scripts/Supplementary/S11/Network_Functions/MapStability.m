if PlotStability==1

    % Choose2NeuronsCA3=[10,11];
    % Parameters
if Context==1
chosen_neurons = Choose2NeuronsCA3; % Replace with indices of chosen neurons
elseif Context==2
chosen_neurons = Choose2NeuronsCA3; % Replace with indices of chosen neurons
end

colors = {'r', 'b'}; % Red for neuron1, Blue for neuron2

% Initialize figure
figure;
hold on;

firing_rate_snapshots_CA3=Rates{Epoch};
% Iterate over chosen neurons
for n = 1:2
    chosen_neuron = chosen_neurons(n); % Current chosen neuron
    color = colors{n}; % Color for the plot

    % Get firing rate snapshots for the chosen neuron
    num_time_points = length(firing_rate_snapshots_CA3);
    correlations = zeros(1, num_time_points); % Store correlations

    % Extract the last time point as reference
    reference_snapshot = firing_rate_snapshots_CA3{end}{chosen_neuron};
    reference_snapshot = reference_snapshot(:); % Flatten into a vector

    % Compute Pearson correlation with each time point
    for t = 1:num_time_points
        current_snapshot = firing_rate_snapshots_CA3{t}{chosen_neuron};
        current_snapshot = current_snapshot(:); % Flatten into a vector
        correlations(t) = corr(reference_snapshot, current_snapshot); % Pearson correlation
    end

    % Plot correlations as a function of time
    plot(1:num_time_points, correlations, 'o', 'Color', color, 'MarkerFaceColor', 'none', 'DisplayName', sprintf('Neuron %d', chosen_neuron));
end

% Customize plot
xlabel('Time Points');
ylabel('Pearson Correlation');
legend('show');
title('Correlation of Last Time Point with All Time Points');
hold off;
ylim([-1,1])
end