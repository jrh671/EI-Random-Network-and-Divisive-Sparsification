

% Load or define your Firing_Rates matrix
Firing_Rates= SlidingAverage_s(TemporalRates{Epoch,1},.1);
% Firing_Rates= TemporalRates{Epoch,1};

% Step 0: Exclude neurons with no activity
% Firing_Rates=Firing_Rates+eps*rand(size(Firing_Rates));
activeNeurons = any(Firing_Rates, 2);
Firing_Rates = Firing_Rates(activeNeurons,:);
% Firing_Rates = zscore(Firing_Rates, 0, 2); % Normalize each row

% Keep track of which neurons are active
activeIndices = find(activeNeurons);
[numCells, numTimePoints] = size(Firing_Rates);

% Original number of cells
originalNumCells = N;

% Step 1: Compute Kendall tau correlation for each cell pair
kendallTauMatrix = corr(Firing_Rates', 'Type', 'Pearson');

% Step 2: Calculate the proportion of significantly negative correlations for each cell
proportionNegativeCorrelations = zeros(numCells, 1);
threshold = -0.05;
clear MeanCorrelations
for i = 1:length(activeIndices)
    negativeCorrelations = kendallTauMatrix(i, :) <= threshold;
    negativeCorrelations(i) = [];  % Remove the self-correlation
    proportionNegativeCorrelations(i) = nansum(negativeCorrelations) / (numCells - 1);
    MeanCorrelations(i) = nanmean(kendallTauMatrix(i, :));

end

% Step 3: Filter MemoryCells to include only those that are still active
% Find indices of activeIndices where the value matches MCells
matchingIndices = find(ismember(activeIndices, MCells));

% Split activeIndices into activeMemoryCells and activeNonMemoryCells based on indices
activeMemoryCells = matchingIndices; % Indices of active memory cells
activeNonMemoryCells = setdiff(1:numel(activeIndices), matchingIndices); % Indices of active non-memory cells

% Step 4: Compute the mean proportion of negative correlations for the active MemoryCells subset
observedMeanProportionMemory = nanmean(proportionNegativeCorrelations(activeMemoryCells));

% Step 5: Compute the mean proportion of negative correlations for active non-MemoryCells subset
observedMeanProportionNonMemory = nanmean(proportionNegativeCorrelations(activeNonMemoryCells));

% Step 6: Generate a distribution of mean proportions for random subsets of the same size
numPermutations = 10000;  % Number of random subsets to generate
randomMeanProportions = zeros(numPermutations, 1);
subsetSize = length(activeMemoryCells);

for i = 1:numPermutations
    randomIndices = randperm(numCells, subsetSize);
    randomMeanProportions(i) = nanmean(proportionNegativeCorrelations(randomIndices));
end

% Step 7: Perform a permutation test to compare Memory vs Non-Memory Cells
observedDifference = observedMeanProportionMemory - observedMeanProportionNonMemory;
permDifference = zeros(numPermutations, 1);

combined = [proportionNegativeCorrelations(activeMemoryCells); proportionNegativeCorrelations(activeNonMemoryCells)];
nMemory = length(activeMemoryCells);
nNonMemory = length(activeNonMemoryCells);

for i = 1:numPermutations
    permIndices = randperm(nMemory + nNonMemory);
    permMemory = combined(permIndices(1:nMemory));
    permNonMemory = combined(permIndices(nMemory+1:end));
    permDifference(i) = nanmean(permMemory) - nanmean(permNonMemory);
end

% Step 1: Sort neurons by their proportion of anti-cofiring
[~, sortedIndices] = sort(proportionNegativeCorrelations, 'descend');

% Step 2: Reorder the Kendall Tau correlation matrix based on sorted indices
sortedKendallTauMatrix = kendallTauMatrix(flipud(sortedIndices), flipud(sortedIndices));

% Step 3: Identify where memory cells appear in the sorted list
sortedMemoryIndices = find(ismember(flipud(sortedIndices), activeMemoryCells));
if Single==0
StoreMemoryCells{Trial}(Context,Epoch)=observedMeanProportionMemory;
StoreNonMemoryCells{Trial}(Context,Epoch)=observedMeanProportionNonMemory;
end

observedDifference

if PlotAntiCo==1
% Step 4: Plot the original correlation matrix with a blue-white-red colormap
figure;
imagesc(sortedKendallTauMatrix);
colormap(redbluecmap); % Custom colormap function for blue-white-red
colorbar;
caxis([-.05, .05]); % Keep the original color scale
hold on;

% Step 5: Overlay a black mask on the upper right triangle
[nRows, nCols] = size(sortedKendallTauMatrix);
[xIdx, yIdx] = find(triu(ones(nRows, nCols), 1)); % Get upper triangle indices

% Use a patch to overlay a fully opaque black mask
for k = 1:length(xIdx)
    x = [yIdx(k)-0.5, yIdx(k)+0.5, yIdx(k)+0.5, yIdx(k)-0.5]; % X-coordinates
    y = [xIdx(k)-0.5, xIdx(k)-0.5, xIdx(k)+0.5, xIdx(k)+0.5]; % Y-coordinates
    fill(x, y, 'k', 'FaceAlpha', 1, 'EdgeColor', 'none'); % Fully opaque black
end

% Step 6: Overlay green lines for memory cell locations
for i = 1:length(sortedMemoryIndices)
    idx = sortedMemoryIndices(i);
    xline(idx, 'g', 'LineWidth', 1.5); % Vertical green line
    yline(idx, 'g', 'LineWidth', 1.5); % Horizontal green line
end

hold off

end


% Custom colormap function for blue-white-red
function cmap = redbluecmap()
    n = 256; % Define colormap resolution
    r = [(0:n/2-1)/(n/2), ones(1,n/2)]'; % Red component
    g = [(0:n/2-1)/(n/2), flip((0:n/2-1)/(n/2))]'; % Green component
    b = [ones(1,n/2), flip((0:n/2-1)/(n/2))]'; % Blue component
    cmap = [r g b]; % Combine into colormap
end


