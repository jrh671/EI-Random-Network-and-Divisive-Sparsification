zeroCells = cellfun(@(x) all(x(:) == 0), M_prim_CA3);

zeroIndices = find(zeroCells);

nonZeroCells = ~zeroCells;         
nonZeroIndices = find(nonZeroCells); % Get the indices of these cells.

chosen_neurons_CA3 = nonZeroIndices;

for W=1:length(nonZeroIndices)
    chosen_neurons_CA3ave = chosen_neurons_CA3(W);
    'Neuron'    
W

    cumulative_activity_map_CA3 = zeros(P, P, length(chosen_neurons_CA3ave));
    occupancy_map = zeros(P, P);


    for t = 1:T

        firing_rates_CA3 = FiringRates{t};


        % Get the current position in the trajectory
        pos = round(trajectory(t, :)); 

        if all(pos > 0) && all(pos <= P) 
            occupancy_map(pos(1), pos(2)) = occupancy_map(pos(1), pos(2)) + 1;

            % Update the cumulative activity map
            for i = 1:N
                if ismember(i, chosen_neurons_CA3ave)
                    neuron_index = find(chosen_neurons_CA3ave == i);
                    cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) = ...
                        cumulative_activity_map_CA3(pos(1), pos(2), neuron_index) + firing_rates_CA3{i}(pos(1), pos(2));
                end
            end

          
        end
    end

    occupancy_map(occupancy_map == 0) = 1; 
    for n = 1
        cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) ./ occupancy_map;
    end

    for n = 1
        max_val = max(cumulative_activity_map_CA3(:, :, n), [], 'all');
        if max_val > 0
            cumulative_activity_map_CA3(:, :, n) = cumulative_activity_map_CA3(:, :, n) / max_val;
        end
    end


matrix = cumulative_activity_map_CA3;
% Find the pixel with the highest activity
[~, maxIndex] = max(matrix(:)); 
[maxRow, maxCol] = ind2sub(size(matrix), maxIndex); 

numRows = size(matrix, 1);
numCols = size(matrix, 2);
maxValue = matrix(maxRow, maxCol); 
rValues = []; 
averageValues = []; 

% Loop through radial distances
for r = 0:max(numRows, numCols)
    % Create a mask for pixels at distance r
    mask = false(size(matrix));
    for row = 1:numRows
        for col = 1:numCols
            % Calculate the radial distance
            distance = sqrt((row - maxRow)^2 + (col - maxCol)^2);
            if abs(distance - r) < 0.5 
                mask(row, col) = true;
            end
        end
    end
    
    valuesAtR = matrix(mask);
    if ~isempty(valuesAtR)
        rValues(end + 1) = r; %#
        averageValues(end + 1) = mean(valuesAtR);
    end
end

xValues = [-fliplr(rValues(2:end)), 0, rValues(2:end)]; 
yValues = [fliplr(averageValues(2:end)), maxValue, averageValues(2:end)]; 


Xall{W}=xValues;
Yall{W}=yValues;
end

