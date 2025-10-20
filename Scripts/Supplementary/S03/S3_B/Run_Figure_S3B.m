addpath('./Data');
addpath('./Code');

%% Choose slow or fast membrane time constant results

TmChoice = 'Slow'; %Slow or Fast Time Constant

%% 1: Plurality | 2: Template | 3: L1| 4: L2| 5: Linear
SaveMean=cell(8,1);
Labels = {'Plurality','Template','L1','L2','Linear'};
Instance=1;
Lam=linspace(0,1,10);
cmap = jet(length(Lam));

Sparsity=0;

figure;
CRange = [0, 1];  

for I = 1:10 %Out of 5

if TmChoice == 'Slow'
load('NoNoise_results_50In_500Ex.mat')
else
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S03/S3_B/Data/resultsFast.mat')
end

Decoder=3;
lambda=Lam(I);
RunDecoder_Stats;
end
title('Small L1 Coefficient (0 to 1)')
clim([0,1])

if Sparsity==1
title('L1 Sparsity');
xlabel('Inhibitory Strength');
ylabel('Sparsity %InActive');
end

Lam=linspace(0,5,10);
cmap = jet(length(Lam));
figure;
CRange = [0, 5];  
for I = 1:10 %Out of 5
Decoder=3;
lambda=Lam(I);
RunDecoder_Stats;
end
title('Large L1 Coefficient (0 to 5)')
clim([0,5])

if Sparsity==1
ylim([0.9,1])
title('L1 Sparsity');
xlabel('Inhibitory Strength');
ylabel('Sparsity %InActive');
end

% SaveMeans is a cell array where each cell is 21x1 double
M = horzcat(SaveMean{:});   % or simply: M = [SaveMeans{:}];
% M is 21 x numel(SaveMean)
