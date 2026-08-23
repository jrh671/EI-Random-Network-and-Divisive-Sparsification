FigureToPlot=1;
PlotFigure=1;
% Convert collected data to matrices
CollectReverseMatrix = cell2mat(CollectReverse');
CollectFreezeMatrix = cell2mat(CollectFreeze');

% Normalize each data set by their respective maximum values
maxReverse = max(CollectReverseMatrix, [], 'all');
maxFreeze = max(CollectFreezeMatrix, [], 'all');

NormalizedReverseMatrix = CollectReverseMatrix;
NormalizedFreezeMatrix = CollectFreezeMatrix;

%% Plot

% figure;
hold on;

% Boxplot 
positions = [1, 2, 4, 5, 7, 8]; 
widths = 0.5; % Set the width of the boxplots 

h1 = boxplot(NormalizedFreezeMatrix, 'Positions', positions(1:2:end), 'Widths', widths, 'Labels', {'F1', 'F2', 'F3'}, 'Color', 'k');
h2 = boxplot(NormalizedReverseMatrix, 'Positions', positions(2:2:end), 'Widths', widths, 'Labels', {'R1', 'R2', 'R3'}, 'Color', 'k');

h1_Box = findobj(h1, 'Tag', 'Box');
h2_Box = findobj(h2, 'Tag', 'Box');
for j = 1:length(h1_Box)
    patch(get(h1_Box(j), 'XData'), get(h1_Box(j), 'YData'), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
end
for j = 1:length(h2_Box)
    patch(get(h2_Box(j), 'XData'), get(h2_Box(j), 'YData'), 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
end

% Overlay actual data points with jitter 
for j = 1:3
    x_jitter_freeze = positions(2*j-1) + 0.15 * (rand(1, numSims) - 0.5);
    scatter(x_jitter_freeze, NormalizedFreezeMatrix(:, j), 'filled', 'k', 'MarkerFaceAlpha', 0.5);

    x_jitter_reverse = positions(2*j) + 0.15 * (rand(1, numSims) - 0.5);
    scatter(x_jitter_reverse, NormalizedReverseMatrix(:, j), 'filled', 'k', 'MarkerFaceAlpha', 0.5);
    
    for idx = 1:numSims
        plot([x_jitter_freeze(idx), x_jitter_reverse(idx)], ...
             [NormalizedFreezeMatrix(idx, j), NormalizedReverseMatrix(idx, j)], ...
             'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]); 
    end
end

title('Paired Boxplots');
ylabel('Normalized Values');
xlabel('Freeze vs Reversal');
xticks(mean(reshape(positions, 2, []), 1)); 
xticklabels({'1', '2', '3'}); 
xlim([0 9]); 
hold off;

if PlotFigure == 1

    FiringRateExa=CollectAvoidRates{1};
    TrainGoal; 

    %% Plot Results
    i=1;

for i = [1 3 4 9 10]
SpecialPos = i; 

figure;
plot(CollectPositions(:, SpecialPos), 'b'); 
hold on; 

yline(SpecialPositions(SpecialPos), 'r--', 'LineWidth', 2); 

title(TaggedSize(i)/n_excit);
xlabel('Index');
ylabel('Value');
hold off;
set(gca, 'YDir', 'reverse');

figure;imagesc(CollectRates{i}(CollectIndices{1,i},:),[0 1]);colorbar;title(TaggedSize(i)/n_excit);


end

    
end
