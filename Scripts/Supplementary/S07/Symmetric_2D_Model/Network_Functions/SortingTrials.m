
if AAB==1

if PreLoadData==1
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S07/Symmetric_2D_Model/PreRun_Data/StoreMemory_AAB.mat')
end

numTrials = numel(StoreMemoryCells);
F_values1 = zeros(1, numTrials);

% Compute F for each trial
for i = 1:numTrials
    trialMatrix1{i} = StoreMemoryCells{i}-StoreNonMemoryCells{i}; % Extract the 2x6 matrix
    F_values1(i) = mean([trialMatrix1{i}(1,4) - trialMatrix1{i}(1,6), trialMatrix1{i}(2,4) - trialMatrix1{i}(2,6)]);
end

% Sort trials by F in descending order
[FSort, sortedIndices_AAB] = sort(F_values1, 'descend');

for i=1:numTrials
TrialData_AAB{i}=trialMatrix1{sortedIndices_AAB(i)};
end


TopPairs_AAB=pairs(sortedIndices_AAB,:);

StoreMemory_AAB = '/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S07/Symmetric_2D_Model/PreRun_Data/StoreMemory_AAB.mat';

save(StoreMemory_AAB, 'TopPairs_AAB', '-append');

else
if PreLoadData==1
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S07/Symmetric_2D_Model/PreRun_Data/StoreMemory_ABA.mat')
end

% Compute F for each trial
for i = 1:numTrials
    trialMatrix2{i} = StoreMemoryCells{i}-StoreNonMemoryCells{i}; % Extract the 2x6 matrix
    F_values2(i) = mean([trialMatrix2{i}(1,4) - trialMatrix2{i}(1,6), trialMatrix2{i}(2,4) - trialMatrix2{i}(2,6)]);
end

% Sort trials by F in descending order
[FSort, sortedIndices_ABA] = sort(F_values2, 'descend');

for i=1:numTrials
TrialData_ABA{i}=trialMatrix2{sortedIndices_ABA(i)};
end

TopPairs_ABA=pairs(sortedIndices_ABA,:);

StoreMemory_ABA = '/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S07/Symmetric_2D_Model/PreRun_Data/StoreMemory_ABA.mat';

save(StoreMemory_ABA, 'TopPairs_ABA', '-append');

end