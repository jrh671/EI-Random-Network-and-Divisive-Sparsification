% Create a logical array where each element is true if the corresponding 10x10 matrix is all zeros.
zeroCells = cellfun(@(x) all(x(:) == 0), M_prim_CA3);

% Optionally, get the indices of matrices that are all zeros:
zeroIndices = find(zeroCells);

nonZeroCells = ~zeroCells;         % Logical array where true indicates the cell is not all zeros.
nonZeroIndices = find(nonZeroCells); % Get the indices of these cells.

chosen_neurons_CA3 = nonZeroIndices;

for W=1:length(nonZeroIndices)
    chosen_neurons_CA3ave = chosen_neurons_CA3(W);
    'Neuron'    
W

    % Initialize arrays to store the cumulative activity at each location for chosen neurons
    cumulative_activity_map_CA3 = zeros(P, P, length(chosen_neurons_CA3ave));
    occupancy_map = zeros(P, P);


    for t = 1:T

        firing_rates_CA3 = FiringRates{t};


        % Get the current position in the trajectory
        pos = round(trajectory(t, :)); % Round to the nearest integer to use as indices

        if all(pos > 0) && all(pos <= P) % Ensure the position is within bounds
            occupancy_map(pos(1), pos(2)) = occupancy_map(pos(1), pos(2)) + 1;

            % Update the cumulative activity map for chosen neurons in CA3
            for i = 1:N
                if ismember(i, chosen_neurons_CA3ave)
                    neuron_index = find(chosen_neurons_CA3ave == i);
                    cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) = ...
                        cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) + firing_rates_CA3{i}(pos(1), pos(2));
                end
            end

          
        end
    end

    % Normalize the activity maps by the occupancy map
    occupancy_map(occupancy_map == 0) = 1; % To avoid division by zero
    for n = 1
        cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) ./ occupancy_map;
    end

    % Normalize the activity maps by the max value of each neuron
    for n = 1
        max_val = max(cumulative_activity_map_CA3(:, :, n), [], 'all');
        if max_val > 0
            cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) / max_val;
        end
    end


% Assuming cumulative_activity_CA3 is a 20x20 matrix
matrix = cumulative_activity_map_CA3;
% Find the pixel with the highest activity
[~, maxIndex] = max(matrix(:)); % Get the linear index of the max value
[maxRow, maxCol] = ind2sub(size(matrix), maxIndex); % Convert to row and column indices

% Initialize variables
numRows = size(matrix, 1);
numCols = size(matrix, 2);
maxValue = matrix(maxRow, maxCol); % Value of the highest activity pixel
rValues = []; % Stores radial distances
averageValues = []; % Stores average pixel values at radius r

% Loop through radial distances
for r = 0:max(numRows, numCols)
    % Create a mask for pixels at distance r
    mask = false(size(matrix));
    for row = 1:numRows
        for col = 1:numCols
            % Calculate the radial distance
            distance = sqrt((row - maxRow)^2 + (col - maxCol)^2);
            if abs(distance - r) < 0.5 % Allowing for some numerical tolerance
                mask(row, col) = true;
            end
        end
    end
    
    % Get the average value of pixels at this distance
    valuesAtR = matrix(mask);
    if ~isempty(valuesAtR)
        rValues(end + 1) = r; %#
        averageValues(end + 1) = mean(valuesAtR);
    end
end

% Prepare data for plotting
xValues = [-fliplr(rValues(2:end)), 0, rValues(2:end)]; % Symmetric radial distances
yValues = [fliplr(averageValues(2:end)), maxValue, averageValues(2:end)]; % Corresponding pixel values


Xall{W}=xValues;
Yall{W}=yValues;
end

