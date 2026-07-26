% Create a logical array where each element is true if the corresponding 10x10 matrix is all zeros.
zeroCells = cellfun(@(x) all(x(:) == 0), M_prim_CA3);

% Optionally, get the indices of matrices that are all zeros:
zeroIndices = find(zeroCells);

nonZeroCells = ~zeroCells;         % Logical array where true indicates the cell is not all zeros.
nonZeroIndices = find(nonZeroCells); % Get the indices of these cells.

length(nonZeroIndices)
RateMap = zeros(P, P, length(nonZeroIndices));

for W=1:length(nonZeroIndices)

    chosen_neurons_CA3ave = nonZeroIndices(W);

    % Initialize arrays to store the cumulative activity at each location for chosen neurons
    occupancy_map = zeros(P, P);

    for t = 1:T

        firing_rates_CA3 = FiringRates{t}{chosen_neurons_CA3ave};

        % Get the current position in the trajectory
        pos = round(trajectory(t, :)); % Round to the nearest integer to use as indices

        if all(pos > 0) && all(pos <= P) % Ensure the position is within bounds
            occupancy_map(pos(1), pos(2)) = occupancy_map(pos(1), pos(2)) + 1;

            % Update the cumulative activity map for chosen neurons in CA3
                    RateMap(pos(1), pos(2), W) = ...
                        RateMap(pos(1), pos(2), W) + firing_rates_CA3(pos(1), pos(2));

          
        end
    end

    % Normalize the activity maps by the occupancy map
    occupancy_map(occupancy_map == 0) = 1; % To avoid division by zero
    for n = 1:length(nonZeroIndices)

        RateMap(:, :, n) = RateMap(:, :, n) ./ occupancy_map;
    end

    % Normalize the activity maps by the max value of each neuron
    for n = 1:length(nonZeroIndices)

        max_val = max(RateMap(:, :, n), [], 'all');
        if max_val > 0
            RateMap(:, :, n) = RateMap(:, :, n) / max_val;
        end
    end

end