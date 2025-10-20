% Assuming FiringRateExb is already defined 
num_neurons = size(FiringRateExb, 1);
num_bins = size(FiringRateExb(:,n_pos+1:end), 2);

% Initialize vectors to store the first time reaching 20 spikes and max spikes
first_reach_20 = nan(num_neurons, 1);
max_spikes = nan(num_neurons, 1);

% Loop through each neuron to find the first time bin reaching 20 spikes and the max spikes
for neuron = 1:num_neurons
    cumulative_spikes = cumsum(FiringRateExb(neuron, n_pos+1:end));
    first_time = find(cumulative_spikes >= 20, 1);
    
    if ~isempty(first_time)
        first_reach_20(neuron, 1) = first_time;
    end
    
    max_spikes(neuron, 1) = max(cumulative_spikes);
end

% Separate indices for neurons that reached 20 spikes and those that didn't
reached_20 = ~isnan(first_reach_20);
not_reached_20 = isnan(first_reach_20);
% Sort neurons that reached 20 spikes by the time they first reached 20 spikes
[~, sort_reached_20] = sort(first_reach_20(reached_20));

% Sort neurons that did not reach 20 spikes by maximum spikes in descending order
[~, sort_not_reached_20] = sort(max_spikes(not_reached_20), 'descend');

% Combine the sorted indices
sorted_indices_reached_20 = find(reached_20);
sorted_indices_not_reached_20 = find(not_reached_20);

% Combine the final sorted indices
Indexical2 = [sorted_indices_reached_20(sort_reached_20); sorted_indices_not_reached_20(sort_not_reached_20)];

