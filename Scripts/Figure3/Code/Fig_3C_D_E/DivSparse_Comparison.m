addpath('./EI_Functions');
addpath('./Helper_Functions');

ShowSparsity=0; %0 For Decoding, 1 for Sparsity
ShowRates=1;

LoadPreRun;

num_thresholds=21;
thresholds = linspace(0, 1, num_thresholds).^(1/10);

DivSparseOn=0;

GenerateActivities;

for Idx = 1:N_params
    Idx
    % Extract 10 samples per Idx
    FiringRateF_samples = X1(Idx, :);  
    FiringRateS_samples = X2(Idx, :);
    Activity_samples = X_DivSparse(Idx, :);

    % Preallocate arrays for sparsity values
    num_samples = length(FiringRateF_samples);  
    SparsityFiringFStore = zeros(num_samples, 1);
    SparsityFiringSStore = zeros(num_samples, 1);
    SparsityActivityStore = zeros(num_samples, 1);

    % Loop through each of the samples
    for sampleIdx = 1:10
        % Extract single sample
        FiringRateF = FiringRateF_samples{sampleIdx};  
        FiringRateS = FiringRateS_samples{sampleIdx};  
        Activity = Activity_samples{sampleIdx};  

        ComputeSparsity;  

        % Store sparsity values
        SparsityFiringFStore(sampleIdx) = SparsityFiringF;  
        SparsityFiringSStore(sampleIdx) = SparsityFiringS;
        SparsityActivityStore(sampleIdx) = SparsityActivity;
    end

    SparsityLIFFmean(Idx,1) = mean(SparsityFiringFStore);
    
    SparsityLIFSmean(Idx,1) = mean(SparsityFiringSStore);
    
    SparsityModmean(Idx,1) = mean(SparsityActivityStore);
end



Overlap=1;
Sparsity=SparsityLIFFmean';
X=X1;
FirstPlot=1;
Run_Figure3_Decoder;
FirstPlot=0;

clear decoding_error

Overlap=2;
Sparsity=SparsityLIFSmean';
X=X2;
Run_Figure3_Decoder;

clear decoding_error

Overlap=3;
DivSparseOn=1;
Sparsity=SparsityModmean';
X=X_DivSparse;

Run_Figure3_Decoder;


numIdx1 = [3, 4, 9]; 
numIdx2 = [1, 11, 19];

if ShowRates==1
for Idx = 1:length(numIdx1)
    

    FiringRateF = X1{numIdx1(Idx), 1};
    FiringRateS = X2{numIdx1(Idx), 1};
    SaveRates{Idx}=FiringRateS(Indices2, 1:end);
    Activity = X_DivSparse{numIdx2(Idx), 1};    
    SaveActivity{Idx}=Activity(neuron_indices, 1:end);



    figure;
    imagesc(FiringRateS(Indices2, 1:end));
    colorbar;
    clim([0, 1]); 
    title(['FiringRate for Idx ', num2str(Idx)]);
    xlabel('Time');
    ylabel('Neuron');


    figure;
    imagesc(Activity(neuron_indices, 1:end));
    colorbar;
    clim([0, 1]); % Adjust the color limits
    title(['Activity for Idx ', num2str(Idx)]);
    xlabel('Time');
    ylabel('Neuron');
end
end