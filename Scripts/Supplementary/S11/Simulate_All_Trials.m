addpath('./Network_Functions');


Single=0;
PlotAntiCo=0;
AAB=0;
GenerateTuning=1;
GenerateActivity=1;
Plasticity=1;
PlotMaps=0;
PlotStability=0;
VisualizePlasticity=0;
VisualizePositions=0;
Rates = cell(8,1);
Sorting = cell(8,1);
MemoryCells = cell(8,1);
TemporalRates = cell(8,1);
TemporalSpikes = cell(8,1);
numContexts=2;
numEpochs=6;

if VisualizePlasticity==1
ContextGoal = [2]; 
else
    ContextGoal = [1,2]; 
end

% ContextGoal = [1,2]; 
% 
% Define range
numbers = 0:10;

% Generate all unique pairs
pairs1 = nchoosek(numbers, 2); % 55x2 matrix
pairs2 = fliplr(nchoosek(numbers, 2)); % 55x2 matrix
pairs=pairs1;%[pairs1',pairs2']';

numTrials=length(pairs); 

% Preallocate cell arrays
StoreMemoryCells = cell(numTrials, 1);
StoreNonMemoryCells = cell(numTrials, 1);

% Initialize each trial with a zero matrix of size (numContexts x numEpochs)
for Trial = 1:numTrials
    StoreMemoryCells{Trial} = zeros(numContexts, numEpochs);
    StoreNonMemoryCells{Trial} = zeros(numContexts, numEpochs);
end

for Trial=1:numTrials

    Trial

for Context = ContextGoal
    clear MCells MemoryCells

    if Context == 1
        Epochs = 1:6;
    elseif Context == 2
        Epochs = 1:6;
    end


    for Epoch = Epochs
    
        % fprintf('Epoch %d\n', Epoch);
    
        % Parameters
        N = 100;         % Number of neurons in CA3
        K = 300;         % Number of inputs
        P = 10;          % Number of positions per dimension (total positions = P x P)
        T = 1200;         % Number of time steps
        StartTime = 50; %End-timee of visualization
        EndTime = 150; %End-time of visualization
        GoalRadius=2;
        
        % Define seed1 and seed2 for each trial
        seed1 = pairs(Trial,1);
        seed2 = pairs(Trial,2);

        noise_std = 0.2; % Uniform Noise
        MemoryLimit=20;

        % Chosen neurons for visualization 
        Choose2NeuronsCA3=[10,48];

        % Define special location and radius
        special_location = [7, 2];
    

            if Epoch == 2 
                special_radius = GoalRadius; 
            else
                special_radius = 0;
            end

        
        %% Run First Without Activity-Dependent-Plasticity
        
        SetParameters4D;
        
         if Context==1
            if Epoch < 3
                if GenerateTuning==1
                        Seed = seed1;           
                    RunGenerateTuning;
                    SaveConnections1 = A_CA3;
                    SaveTuning1 = B;
                end
            elseif Epoch == 3 || Epoch == 4
                if AAB==1
                    Seed = seed1;
                else
                    Seed = seed2;
                end
                    Remap_InputTuning;
                    SaveConnections2 = A_CA3;
                    SaveTuning2 = B;
            elseif Epoch == 5 || Epoch == 6
                    if AAB==1
                        Seed = seed2;
                    else
                        Seed = seed1;
                    end
                Remap_InputTuning;
                    SaveTuning3 = B;
            end

        elseif Context==2
            if Epoch < 3
                if GenerateTuning==1
                    Seed = seed2;
                    RunGenerateTuning;
                    SaveConnections1 = A_CA3;
                    SaveTuning1 = B;
                end
            elseif Epoch == 3 || Epoch == 4
                    if AAB==1
                        Seed = seed2;
                    else
                        Seed = seed1;
                    end
                Remap_InputTuning;
                    SaveConnections2 = A_CA3;
                    SaveTuning2 = B;
            elseif Epoch == 5 || Epoch == 6
                    if AAB==1
                        Seed = seed1;
                    else
                        Seed = seed2;
                    end
                Remap_InputTuning;
                    SaveTuning3 = B;
            end
        end

        
        if GenerateActivity==1
            RunGenerateActivity;
            Rates{Epoch} = firing_rate_snapshots_CA3;
            Sorting{Epoch} = sorted_indices_CA3;
            MemoryCells{Epoch} = memory_neurons_CA3;
        end
    
    end

    Visualize2DMaps;

end
end

PreLoadData=0;
MemoryStoreAntiCofiring_Diff
