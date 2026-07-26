addpath('./Data');
addpath('./Code');

%% Choose Slow Or Fast Time Constant Results

TmChoice='Fast'; %Fast or Slow

%% Colors PreRun
% Colors = ['r','b','k','g','c']; IDs = [Assembly Tagging, Template Matchiing, Lasso, Ridge, Unregularized Linear]

%% Labels 1: Plurality | 2: Template | 3: L1| 4: L2| 5: Linear 

Labels = {'Plurality','Template','L1','L2','Linear'};
figure
SaveMean=cell(8,1);

for Decoder = 1:5 %Out of 5

if strcmp(TmChoice,'Slow')

    % Sparse/full range: 0 to 0.3, 21 points
    S1 = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsSlow1.mat');


    % Dense short range: 0 to 0.1
    S2 = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsSlow2.mat')

    results_sparse = S1.results;
    results_dense  = S2.results;

    conditions_sparse = linspace(0,0.3,21);
    conditions_dense  = linspace(0,0.1,21);

    % Keep all dense points from 0 to 0.1
    keep_dense = conditions_dense <= 0.1;

    % Keep sparse points after 0.1 only, avoiding duplicate 0.1
    keep_sparse = conditions_sparse > 0.1;

    combined_conditions = [conditions_dense(keep_dense), conditions_sparse(keep_sparse)];

    % Each condition has 10 samples
    nSamples = 10;

    combined_results = {};

    counter = 1;

    for n = find(keep_dense)
        for s = 1:nSamples
            old_idx = (n-1)*nSamples + s;
            combined_results{counter} = results_dense{old_idx};
            counter = counter + 1;
        end
    end

    for n = find(keep_sparse)
        for s = 1:nSamples
            old_idx = (n-1)*nSamples + s;
            combined_results{counter} = results_sparse{old_idx};
            counter = counter + 1;
        end
    end

    results = combined_results;
    conditions = combined_conditions;

else


    % Fast sparse range: 0 to 0.5, step 0.025
    F_sparse = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsFast1.mat');
    
    % Fast dense range: 0 to 0.3, step 0.015
    F_dense = load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2C_2D/2C/PreRun_Data/resultsFast2.mat');

    results_dense  = F_dense.results;
    results_sparse = F_sparse.results;

    nSamples = 10;

    OhmVals = 1:21;

    conditions_dense  = (OhmVals - 1) * 0.015;  % 0 to 0.3
    conditions_sparse = (OhmVals - 1) * 0.025;  % 0 to 0.5

    % Keep dense through 0.3
    keep_dense = conditions_dense <= 0.3 + eps;

    % Keep sparse only AFTER 0.3, no duplicate 0.3
    keep_sparse = conditions_sparse > 0.3 + eps;

    combined_conditions = [conditions_dense(keep_dense), conditions_sparse(keep_sparse)];

    combined_results = {};
    counter = 1;

    for n = find(keep_dense)
        for s = 1:nSamples
            old_idx = (n-1)*nSamples + s;
            combined_results{counter} = results_dense{old_idx};
            counter = counter + 1;
        end
    end

    for n = find(keep_sparse)
        for s = 1:nSamples
            old_idx = (n-1)*nSamples + s;
            combined_results{counter} = results_sparse{old_idx};
            counter = counter + 1;
        end
    end

    results = combined_results;
    conditions = combined_conditions;

end


Instance=1;
RunDecoder_Stats;
end

% legend(hErrorbar, Labels); % Add legend after all plots are created
