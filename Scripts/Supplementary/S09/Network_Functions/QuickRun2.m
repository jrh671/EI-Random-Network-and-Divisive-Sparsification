
if Average==1
% Calculate the mean and SEM
numNeurons = length(Yall);

% Assuming all X1 arrays are identical (shared radial distances)
[~, maxIdx] = max(cellfun(@length, Xall)); % Find index of longest cell element
X1 = Xall{maxIdx}; % Extract the longest cell element

% Determine the maximum length of vectors in Y1all
maxLength = nanmax(cellfun(@length, Yall));

% Pad each Y1 vector with NaN to make all vectors the same length
Y1Padded = cellfun(@(y) [y, nan(1, maxLength - length(y))], Yall, 'UniformOutput', false);

% Convert the padded cell array to a matrix
Y1Matrix = cell2mat(Y1Padded'); % Each column corresponds to a neuron

% Calculate the mean and SEM across neurons
Y1Mean = nanmean(Y1Matrix, 1);
Y1SEM = nanstd(Y1Matrix, 0, 1) / sqrt(numNeurons);

end

%% Run Next With Activity-Dependent-Plasticity
Plasticity=1;

SetParameters4D;
if GenerateTuning==1
    RunGenerateTuning;
end

if GenerateActivity==1
    'Begin Activity'
    RunGenerateActivity;
end

FiringRates=FiringRates2;
Traj2=trajectory;
FRAve2=cumulative_activity_map_CA3;
AllCells_RateMap;
RateMap2=RateMap;

% Visualization
Visualize2DMaps;



CalcOverdispersion;

if Average == 1
% Calculate the mean and SEM
numNeurons = length(Yall);

% Assuming all X1 arrays are identical (shared radial distances)
[~, maxIdx] = max(cellfun(@length, Xall)); % Find index of longest cell element
X2 = Xall{maxIdx}; % Extract the longest cell element

% Determine the maximum length of vectors in Y1all
maxLength = nanmax(cellfun(@length, Yall));

% Pad each Y1 vector with NaN to make all vectors the same length
Y2Padded = cellfun(@(y) [y, nan(1, maxLength - length(y))], Yall, 'UniformOutput', false);

% Convert the padded cell array to a matrix
Y2Matrix = cell2mat(Y2Padded'); % Each column corresponds to a neuron

% Calculate the mean and SEM across neurons
Y2Mean = nanmean(Y2Matrix, 1);
Y2SEM = nanstd(Y2Matrix, 0, 1) / sqrt(numNeurons);

% Plot the mean with SEM as error bars
figure;
errorbar(X1, Y1Mean(1:length(X1)), Y1SEM(1:length(X1)), 'b-o', 'LineWidth', 1.5);
hold on
errorbar(X2, Y2Mean(1:length(X2)), Y2SEM(1:length(X2)), 'r-o', 'LineWidth', 1.5);
xlabel('Radial Distance (r)');
ylabel('Pixel Value');
title('Mean Pixel Value with SEM Across Neurons');

SheetPFSpreadMean1=[X1',Y1Mean(1:length(X1))'];
SheetPFSpreadSTD1=[X1',Y1SEM(1:length(X1))'];

SheetPFSpreadMean2=[X1',Y2Mean(1:length(X1))'];
SheetPFSpreadSTD2=[X1',Y2SEM(1:length(X1))'];

end
