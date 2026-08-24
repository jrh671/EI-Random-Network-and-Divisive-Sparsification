

Firing_Rates= SlidingAverage_s(TemporalRates{Epoch,1},.1);


activeNeurons = any(Firing_Rates, 2);
Firing_Rates = Firing_Rates(activeNeurons,:);

activeIndices = find(activeNeurons);
[numCells, numTimePoints] = size(Firing_Rates);

originalNumCells = N;

% Compute Correlations
PearsonMatrix = corr(Firing_Rates', 'Type', 'Pearson');

proportionNegativeCorrelations = zeros(numCells, 1);
threshold = -0.05;
clear MeanCorrelations
for i = 1:length(activeIndices)
    negativeCorrelations = PearsonMatrix(i, :) <= threshold;
    negativeCorrelations(i) = [];  
    proportionNegativeCorrelations(i) = nansum(negativeCorrelations) / (numCells - 1);
    MeanCorrelations(i) = nanmean(PearsonMatrix(i, :));

end

matchingIndices = find(ismember(activeIndices, MCells));

% Split activeIndices into activeMemoryCells and activeNonMemoryCells based on indices
activeMemoryCells = matchingIndices; % Indices of active memory cells
activeNonMemoryCells = setdiff(1:numel(activeIndices), matchingIndices); % Indices of active non-memory cells

observedMeanProportionMemory = nanmean(proportionNegativeCorrelations(activeMemoryCells));

observedMeanProportionNonMemory = nanmean(proportionNegativeCorrelations(activeNonMemoryCells));

% Generate a distribution of mean proportions for random subsets
numPermutations = 10000;  % Number of random subsets to generate
randomMeanProportions = zeros(numPermutations, 1);
subsetSize = length(activeMemoryCells);

for i = 1:numPermutations
    randomIndices = randperm(numCells, subsetSize);
    randomMeanProportions(i) = nanmean(proportionNegativeCorrelations(randomIndices));
end

% Permutation test to compare Memory vs Non-Memory Cells
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

[~, sortedIndices] = sort(proportionNegativeCorrelations, 'descend');

sortedPearsonMatrix = PearsonMatrix(flipud(sortedIndices), flipud(sortedIndices));

sortedMemoryIndices = find(ismember(flipud(sortedIndices), activeMemoryCells));
if Single==0
StoreMemoryCells{Trial}(Context,Epoch)=observedMeanProportionMemory;
StoreNonMemoryCells{Trial}(Context,Epoch)=observedMeanProportionNonMemory;
end

observedDifference

if PlotAntiCo==1
figure;
imagesc(sortedPearsonMatrix);
colormap(redbluecmap); 
colorbar;
caxis([-.05, .05]); 
hold on;

[nRows, nCols] = size(sortedPearsonMatrix);
[xIdx, yIdx] = find(triu(ones(nRows, nCols), 1)); 

for k = 1:length(xIdx)
    x = [yIdx(k)-0.5, yIdx(k)+0.5, yIdx(k)+0.5, yIdx(k)-0.5]; % X-coordinates
    y = [xIdx(k)-0.5, xIdx(k)-0.5, xIdx(k)+0.5, xIdx(k)+0.5]; % Y-coordinates
    fill(x, y, 'k', 'FaceAlpha', 1, 'EdgeColor', 'none'); % Fully opaque black
end

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


