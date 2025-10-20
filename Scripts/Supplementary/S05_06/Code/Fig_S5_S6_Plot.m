addpath('Results/')
load('./Results/Results.mat')
addpath("EI_Functions/")
addpath("Helper_Functions/")

%Plots Results

PlotFigure = 1;
FigureToPlot = 9; % Simulation 1 and 9 are in Figure S6 

% Specify the indices of the files you want to use (choose any 10 out of 12)
selectedIndices = randperm(10,10);  % Choose 10 Separate Simulations

% Initialize cell arrays to collect data
numFiles = length(selectedIndices);

% Convert collected data to matrices
CollectReverseMatrix = cell2mat(CollectReverse');
CollectFreezeMatrix = cell2mat(CollectFreeze');

% Normalize each data set by their respective maximum values
maxReverse = max(CollectReverseMatrix, [], 'all');
maxFreeze = max(CollectFreezeMatrix, [], 'all');

NormalizedReverseMatrix = CollectReverseMatrix;% / maxReverse;
NormalizedFreezeMatrix = CollectFreezeMatrix;% / maxFreeze;

% Create a single figure for paired boxplots
figure;
hold on;

% Boxplot positions for grouping (switching order within pairs)
positions = [1, 2, 4, 5, 7, 8]; % Adjust positions for grouping
widths = 0.5; % Set the width of the boxplots to be wider

% Plot boxplots for NormalizedFreezeMatrix first and then NormalizedReverseMatrix
h1 = boxplot(NormalizedFreezeMatrix, 'Positions', positions(1:2:end), 'Widths', widths, 'Labels', {'F1', 'F2', 'F3'}, 'Color', 'k');
h2 = boxplot(NormalizedReverseMatrix, 'Positions', positions(2:2:end), 'Widths', widths, 'Labels', {'R1', 'R2', 'R3'}, 'Color', 'k');

% Set the transparency for the boxes
h1_Box = findobj(h1, 'Tag', 'Box');
h2_Box = findobj(h2, 'Tag', 'Box');
for j = 1:length(h1_Box)
    patch(get(h1_Box(j), 'XData'), get(h1_Box(j), 'YData'), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Red transparent box
end
for j = 1:length(h2_Box)
    patch(get(h2_Box(j), 'XData'), get(h2_Box(j), 'YData'), 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Blue transparent box
end



% Overlay actual data points with jitter and connect them
for j = 1:3
    % Jitter for NormalizedFreezeMatrix
    x_jitter_freeze = positions(2*j-1) + 0.15 * (rand(1, numFiles) - 0.5);
    scatter(x_jitter_freeze, NormalizedFreezeMatrix(:, j), 'filled', 'k', 'MarkerFaceAlpha', 0.5);

    % Jitter for NormalizedReverseMatrix
    x_jitter_reverse = positions(2*j) + 0.15 * (rand(1, numFiles) - 0.5);
    scatter(x_jitter_reverse, NormalizedReverseMatrix(:, j), 'filled', 'k', 'MarkerFaceAlpha', 0.5);
    
    % Connect the data points with a line
    for idx = 1:numFiles
        plot([x_jitter_freeze(idx), x_jitter_reverse(idx)], ...
             [NormalizedFreezeMatrix(idx, j), NormalizedReverseMatrix(idx, j)], ...
             'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]); % Light gray line with transparency
    end
end

SheetBoxFreeze=NormalizedFreezeMatrix;
SheetBoxReverse=NormalizedReverseMatrix;

title('Paired Boxplots');
ylabel('Normalized Values');
xlabel('Components');
xticks(mean(reshape(positions, 2, []), 1)); % Set xticks to the center of pairs
xticklabels({'1', '2', '3'}); % Set xticklabels to match component pairs
xlim([0 9]); % Adjust x limits for better visualization
hold off;

% Adjust figure
sgtitle('Boxplots of CollectReverse and CollectFreeze');

if PlotFigure == 1
    
    F=FigureToPlot;
    FiringRateExa=CollectRates{F};
    neuron_indices=CollectIndices{1,F};

    figure;imagesc(FiringRateExa(neuron_indices,:),[0 1]);colorbar
    SheetActivities=FiringRateExa(neuron_indices,:);

    PosRecall = CollectPositions(:, F);
    figure; plot(PosRecall); ylim([0, 30])
    SheetPositions=PosRecall;
    
    Test_Goals




    Test_Primacy;
end
