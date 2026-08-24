if PlotStability==1

    % Parameters
if Context==1
chosen_neurons = Choose2NeuronsCA3; 
elseif Context==2
chosen_neurons = Choose2NeuronsCA3; 
end

colors = {'r', 'b'}; 
figure;
hold on;

firing_rate_snapshots_CA3=Rates{Epoch};
% Iterate over chosen neurons
for n = 1:2
    chosen_neuron = chosen_neurons(n); % Current chosen neuron
    color = colors{n}; 

    % Get firing rate snapshots 
    num_time_points = length(firing_rate_snapshots_CA3);
    correlations = zeros(1, num_time_points); 

    % Extract the last time point as reference
    reference_snapshot = firing_rate_snapshots_CA3{end}{chosen_neuron};
    reference_snapshot = reference_snapshot(:); 

    % Compute Pearson correlation with each time point
    for t = 1:num_time_points
        current_snapshot = firing_rate_snapshots_CA3{t}{chosen_neuron};
        current_snapshot = current_snapshot(:); 
        correlations(t) = corr(reference_snapshot, current_snapshot);
    end

    plot(1:num_time_points, correlations, 'o', 'Color', color, 'MarkerFaceColor', 'none', 'DisplayName', sprintf('Neuron %d', chosen_neuron));
    SaveCorrelations{n}=correlations';
end

% Customize plot
xlabel('Time Points');
ylabel('Pearson Correlation');
legend('show');
title('Correlation of Last Time Point with All Time Points');
hold off;
ylim([-1,1])
end