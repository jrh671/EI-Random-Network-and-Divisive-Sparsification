addpath('./Data');
addpath('./Code');

%% Choose slow or fast membrane time constant results

TmChoice = 'Slow'; % Slow or Fast Time Constant

%% 1: Plurality | 2: Template | 3: L1 | 4: L2 | 5: Linear

SaveMean = cell(8,1);
Labels = {'Plurality','Template','L1','L2','Linear'};
Instance = 1;
Lam = linspace(0,1,10);
cmap = jet(length(Lam));

Sparsity = 0;

if strcmp(TmChoice,'Slow')

    S1 = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsSlow1.mat');
    S2 = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsSlow2.mat');

    results_sparse = S1.results;
    results_dense  = S2.results;

    OhmVals = 1:21;
    conditions_sparse = (OhmVals - 1) * 0.015;  % 0 to 0.3
    conditions_dense  = (OhmVals - 1) * 0.005;  % 0 to 0.1

    cutoff = 0.1;

elseif strcmp(TmChoice,'Fast')

    F1 = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsFast1.mat');
    F2 = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsFast2.mat');

    results_dense  = F1.results;   % 0 to 0.3, step 0.015
    results_sparse = F2.results;   % 0 to 0.5, step 0.025

    OhmVals = 1:21;
    conditions_dense  = (OhmVals - 1) * 0.015;  % 0 to 0.3
    conditions_sparse = (OhmVals - 1) * 0.025;  % 0 to 0.5

    cutoff = 0.3;

end

keep_dense = conditions_dense <= cutoff + eps;
keep_sparse = conditions_sparse > cutoff + eps;

conditions = [conditions_dense(keep_dense), conditions_sparse(keep_sparse)];

nSamples = 10;
results = {};
counter = 1;

for n = find(keep_dense)
    for s = 1:nSamples
        old_idx = (n-1)*nSamples + s;
        results{counter} = results_dense{old_idx};
        counter = counter + 1;
    end
end

for n = find(keep_sparse)
    for s = 1:nSamples
        old_idx = (n-1)*nSamples + s;
        results{counter} = results_sparse{old_idx};
        counter = counter + 1;
    end
end

figure;
CRange = [0, 1];

for I = 1:10

    Decoder = 3;
    lambda = Lam(I);
    RunDecoder_Stats;

end

title('Small L1 Coefficient (0 to 1)');
clim([0,1]);

if Sparsity == 1
    title('L1 Sparsity');
    xlabel('Inhibitory Strength');
    ylabel('Sparsity %InActive');
end

Lam = linspace(0,5,10);
cmap = jet(length(Lam));

figure;
CRange = [0, 5];

for I = 1:10

    Decoder = 3;
    lambda = Lam(I);
    RunDecoder_Stats;

end

title('Large L1 Coefficient (1 to 5)');
clim([1,5]);

if Sparsity == 1
    ylim([0.9,1]);
    title('L1 Sparsity');
    xlabel('Inhibitory Strength');
    ylabel('Sparsity %InActive');
end

M = horzcat(SaveMean{:});

SaveMeanM = cell2mat(SaveMean');
SaveSTDM = cell2mat(SaveError');
