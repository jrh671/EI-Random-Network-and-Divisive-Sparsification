addpath('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2A_2B/Code/Additional_Simulations/');
addpath('EI_Functions');
addpath('Helper_Functions');




PlotFigure = 1;

% Specify the indices of the simulations
selectedIndices = randperm(10,10);  % Choose 10 Separate Simulations
FigureToPlot = 10;

selectedIndices

numFiles = length(selectedIndices);
CollectReverse = cell(1, numFiles);
CollectFreeze = cell(1, numFiles);
CollectPositions = [];

% Loop over the selected files
for idx = 1:numFiles
    i = selectedIndices(idx);  % Get the actual index from selectedIndices

    % Create the file path for each .mat file
    filePath = sprintf('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2A_2B/Code/Additional_Simulations/%d.mat', i);

    % Load the data
    load(filePath, 'ReverseData', 'FreezeData', 'PosRecall');
    ReverseData(ReverseData < 0) = 0;
    FreezeData(FreezeData < 0) = 0;

    % Collect data
    CollectReverse{idx} = ReverseData;
    CollectFreeze{idx} = FreezeData;
    CollectPositions(:, idx) = PosRecall;
end

CollectReverseMatrix = cell2mat(CollectReverse');
CollectFreezeMatrix = cell2mat(CollectFreeze');

maxReverse = max(CollectReverseMatrix, [], 'all');
maxFreeze = max(CollectFreezeMatrix, [], 'all');

NormalizedReverseMatrix = CollectReverseMatrix / maxReverse;
NormalizedFreezeMatrix = CollectFreezeMatrix / maxFreeze;

%% Plot
figure;
hold on;

% Boxplot 
positions = [1, 2, 4, 5, 7, 8]; 
widths = 0.5; % Set the width of the boxplots

h1 = boxplot(NormalizedFreezeMatrix, 'Positions', positions(1:2:end), 'Widths', widths, 'Labels', {'F1', 'F2', 'F3'}, 'Color', 'k');
h2 = boxplot(NormalizedReverseMatrix, 'Positions', positions(2:2:end), 'Widths', widths, 'Labels', {'R1', 'R2', 'R3'}, 'Color', 'k');

h1_Box = findobj(h1, 'Tag', 'Box');
h2_Box = findobj(h2, 'Tag', 'Box');
for j = 1:length(h1_Box)
    patch(get(h1_Box(j), 'XData'), get(h1_Box(j), 'YData'), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Red transparent box
end
for j = 1:length(h2_Box)
    patch(get(h2_Box(j), 'XData'), get(h2_Box(j), 'YData'), 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Blue transparent box
end

% Overlay actual data points with jitter
for j = 1:3
    x_jitter_freeze = positions(2*j-1) + 0.15 * (rand(1, numFiles) - 0.5);
    scatter(x_jitter_freeze, NormalizedFreezeMatrix(:, j), 'filled', 'k', 'MarkerFaceAlpha', 0.5);

    x_jitter_reverse = positions(2*j) + 0.15 * (rand(1, numFiles) - 0.5);
    scatter(x_jitter_reverse, NormalizedReverseMatrix(:, j), 'filled', 'k', 'MarkerFaceAlpha', 0.5);
    
    for idx = 1:numFiles
        plot([x_jitter_freeze(idx), x_jitter_reverse(idx)], ...
             [NormalizedFreezeMatrix(idx, j), NormalizedReverseMatrix(idx, j)], ...
             'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]); 
    end
end

title('Paired Boxplots of Normalized CollectReverse and CollectFreeze');
ylabel('Normalized Values');
xlabel('Components');
xticks(mean(reshape(positions, 2, []), 1)); 
xticklabels({'1', '2', '3'}); 
xlim([0 9]); 
hold off;

sgtitle('Boxplots of CollectReverse and CollectFreeze');

if PlotFigure == 1
    CollectReverse = cell(1, 10);
    CollectFreeze = cell(1, 10);

    % Loop over the files
    for i = FigureToPlot

       % Create the file path for each .mat file
       filePath = sprintf('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2A_2B/Code/Additional_Simulations/%d.mat', i);
        
        % Load the data
        load(filePath);
        ReverseData(ReverseData < 0) = 0;
        FreezeData(FreezeData < 0) = 0;

        % Collect data
        CollectReverse{i} = ReverseData;
        CollectFreeze{i} = FreezeData;
        CollectPositions(:, i) = PosRecall;
    end

    TrainGoal; 

    %% Plot Results

    PosRecall = CollectPositions(:, i);
    figure; plot(PosRecall); ylim([0, 30])
    pf_cell=pf_cell1;
    Run_EffectiveTuning;
    figure;imagesc(FR(neuron_indices,:),[0 1]);colorbar
    
end
