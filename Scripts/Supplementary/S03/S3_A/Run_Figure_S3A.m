addpath('./Data');
addpath('./Code');

%% Choose Slow Or Fast Time Constant Results

TmChoice='Slow'; %Fast or Slow

%% Colors PreRun
% Colors = ['r','b','k','g','c']; IDs = [Assembly Tagging, Template Matchiing, Lasso, Ridge, Unregularized Linear]

%% Labels 1: Plurality | 2: Template | 3: L1| 4: L2| 5: Linear 

Labels = {'Plurality','Template','L1','L2','Linear'};
figure
SaveMean=cell(8,1);

for Decoder = 1:5 %Out of 5

if TmChoice == 'Slow'
load('NoNoise_results_50In_500Ex.mat')
else
load('resultsFast.mat')
end

Instance=1;
RunDecoder_Stats;
end

% legend(hErrorbar, Labels); % Add legend after all plots are created
