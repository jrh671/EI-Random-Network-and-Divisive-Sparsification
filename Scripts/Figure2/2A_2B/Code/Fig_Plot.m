load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2A_2B/Code/PreRun_Data/Results.mat')
addpath("EI_Functions/")
addpath("Helper_Functions/")

%Plots Results

PlotFigure = 1;
FigureToPlot = 9; % Simulation 1 and 9 are in Figure S6 


selectedIndices = randperm(10,10);  % Choose 10 Separate Simulations
numFiles = length(selectedIndices);

% Convert collected data to matrices
CollectReverseMatrix = cell2mat(CollectReverse');
CollectFreezeMatrix = cell2mat(CollectFreeze');

maxReverse = max(CollectReverseMatrix, [], 'all');
maxFreeze = max(CollectFreezeMatrix, [], 'all');

NormalizedReverseMatrix = CollectReverseMatrix;% / maxReverse;
NormalizedFreezeMatrix = CollectFreezeMatrix;% / maxFreeze;

%% Plot
figure;
hold on;

% Boxplot 
positions = [1, 2, 4, 5, 7, 8]; 
widths = 0.5; % Set the width of the boxplots

% Plot boxplots 
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
             'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]); % Light gray line with transparency
    end
end

SheetBoxFreeze=NormalizedFreezeMatrix;
SheetBoxReverse=NormalizedReverseMatrix;

title('Paired Boxplots');
ylabel('Normalized Values');
xlabel('Components');
xticks(mean(reshape(positions, 2, []), 1)); 
xticklabels({'1', '2', '3'}); 
xlim([0 9]); 
hold off;

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

    load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2A_2B/Code/Additional_Simulations/1.mat')

    Test_Goals




    Test_Primacy;
end
