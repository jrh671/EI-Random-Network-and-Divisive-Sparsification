addpath('./EI_Functions');
addpath('./Helper_Functions');

%% Be Sure To Edit Directory In "LoadPreRun" in Helper_Functions Folder.

ShowSparsity=0; %0 FOr Decoding, 1 for Sparsity
ShowRates=1;

LoadPreRun;

num_thresholds=21;
thresholds = linspace(0, 1, num_thresholds).^(1/10);

NormThresh=0;

GenerateActivities;

for Idx = 1:N_params
    Idx
    % Extract 10 samples per Idx
    FiringRateF_samples = X1(Idx, :);  % Cell array of 10 samples
    FiringRateS_samples = X2(Idx, :);
    Activity_samples = X_NormThresh(Idx, :);

    % Preallocate arrays for sparsity values
    num_samples = length(FiringRateF_samples);  % Should be 10
    SparsityFiringFStore = zeros(num_samples, 1);
    SparsityFiringSStore = zeros(num_samples, 1);
    SparsityActivityStore = zeros(num_samples, 1);

    % Loop through each of the 10 samples
    for sampleIdx = 1:10
        % Extract single sample
        FiringRateF = FiringRateF_samples{sampleIdx};  % Single sample
        FiringRateS = FiringRateS_samples{sampleIdx};  % Single sample
        Activity = Activity_samples{sampleIdx};  % Single sample

        % Compute sparsity for this sample
        ComputeSparsity;  

        % Store sparsity values
        SparsityFiringFStore(sampleIdx) = SparsityFiringF;  % Store computed sparsity
        SparsityFiringSStore(sampleIdx) = SparsityFiringS;
        SparsityActivityStore(sampleIdx) = SparsityActivity;
    end

    % Compute mean for each Idx
    SparsityLIFFmean(Idx,1) = mean(SparsityFiringFStore);
    
    SparsityLIFSmean(Idx,1) = mean(SparsityFiringSStore);
    
    SparsityModmean(Idx,1) = mean(SparsityActivityStore);
end



Overlap=1;
Sparsity=SparsityLIFFmean';
X=X1;
FirstPlot=1;
Run_Figure2_Decoder;
FirstPlot=0;

clear decoding_error

Overlap=2;
Sparsity=SparsityLIFSmean';
X=X2;
Run_Figure2_Decoder;

clear decoding_error

Overlap=3;
NormThresh=1;
Sparsity=SparsityModmean';
X=X_NormThresh;

Run_Figure2_Decoder;


numIdx1 = [3, 4, 9]; % Total number of Idx
numIdx2 = [1, 11, 19]; % Total number of Idx

if ShowRates==1
for Idx = 1:length(numIdx1)
    % Extract FiringRate and Activity for the current Idx
    FiringRateF = X1{numIdx1(Idx), 1};
    FiringRateS = X2{numIdx1(Idx), 1};
    SaveRates{Idx}=FiringRateS(Indices2, 1:end);
    Activity = X_NormThresh{numIdx2(Idx), 1};    
    SaveActivity{Idx}=Activity(neuron_indices, 1:end);



    % Plot FiringRate
    figure;
    imagesc(FiringRateS(Indices2, 1:end));
    colorbar;
    clim([0, 1]); % Adjust the color limits
    title(['FiringRate for Idx ', num2str(Idx)]);
    xlabel('Time');
    ylabel('Neuron');


    % Plot Activity
    figure;
    imagesc(Activity(neuron_indices, 1:end));
    colorbar;
    clim([0, 1]); % Adjust the color limits
    title(['Activity for Idx ', num2str(Idx)]);
    xlabel('Time');
    ylabel('Neuron');
end
end