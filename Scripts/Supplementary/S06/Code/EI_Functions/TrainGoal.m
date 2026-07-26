%Data, ignoring first lap 
FF = FiringRateExa(:, n_pos:end);

data = FF;  

% Define the threshold value for significant peaks
threshold = 22;  

% Define the window size (in number of columns) to detect peaks
windowSize = 30;  

% Initialize a logical matrix to track neurons with peaks within the window
peakMatrix = false(size(data));

% Loop through each neuron (row)
for i = 1:size(data, 1)
    % Loop through each time point with the sliding window
    for j = 1:(size(data, 2) - windowSize + 1)
        % Extract the current window segment
        windowSegment = data(i, j:j + windowSize - 1);
        
        % Check if any value in the window segment exceeds the threshold
        if any(windowSegment > threshold)
            peakMatrix(i, j:j + windowSize - 1) = true;
        end
    end
end

% Filter neurons based on their activity peaks
filteredNeurons = any(peakMatrix, 2);

% Z-score normalization
dataZ = zscore(data, 0, 2);  % Z-score each row (neuron) independently

% Separate filtered and non-chosen neurons
filteredNeuronDataZ = dataZ(filteredNeurons, :);
nonFilteredNeuronDataZ = dataZ(~filteredNeurons, :);

% Calculate the average activity for filtered and non-chosen neurons
averageFiltered = mean(filteredNeuronDataZ, 1);
averageNonFiltered = mean(nonFilteredNeuronDataZ, 1);

if J==1
% Plot the results
figure;
hold on;
plot(averageFiltered, 'b', 'LineWidth', 2);
plot(averageNonFiltered, 'r', 'LineWidth', 2);
title('Average Activity of Filtered and Non-Filtered Neurons');
xlabel('Time Points');
ylabel('Z-Scored Activity');
legend({'Filtered Neurons', 'Non-Filtered Neurons'});
hold off;

Sheet3a=[averageFiltered(:),averageNonFiltered(:)];

end

NumTagged=size(filteredNeuronDataZ,1);


%Collect the tagged neurons
TopNeurons=find(filteredNeurons==1);