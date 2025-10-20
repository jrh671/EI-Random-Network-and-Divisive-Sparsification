addpath('./Network_Functions');
addpath('./PreRun_Data');

AAB=1; %1 = No intermediate remapping. 0 = With intermediate remapping
PairID=54; %% 1 For Fig S9a. #54 for Fig S9d.

if AAB==1
    load('./PreRun_Data/StoreMemory_AAB.mat', 'TopPairs_AAB'); 
    TopPairs=TopPairs_AAB;
else
    load('./PreRun_Data/StoreMemory_ABA.mat', 'TopPairs_ABA'); 
    TopPairs=TopPairs_ABA;
end

Single=1;
PlotAntiCo=1;
PlotStability=1;

GenerateTuning=1;
GenerateActivity=1;
Plasticity=1;
PlotMaps=0;
VisualizePlasticity=0;
VisualizePositions=0;
Rates = cell(8,1);
Sorting = cell(8,1);
MemoryCells = cell(8,1);
TemporalRates = cell(8,1);
TemporalSpikes = cell(8,1);


if VisualizePlasticity==1
    ContextGoal = [1]; 
else
    ContextGoal = [1,2]; 
end




for Context = ContextGoal
    clear MCells MemoryCells

    if Context == 1
        Epochs = 1:6;
    elseif Context == 2
        Epochs = 1:6;
    end


    for Epoch = Epochs
    
        fprintf('Epoch %d\n', Epoch);
    
        % Parameters
        N = 100;         % Number of neurons in CA3
        K = 300;         % Number of inputs
        P = 10;          % Number of positions per dimension (total positions = P x P)
        T = 1200;         % Number of time steps
        StartTime = 10; %End-time of visualization
        EndTime = 100; %End-time of visualization
        GoalRadius=2;
        seed1=TopPairs(PairID,1); %Environment 0
        seed2=TopPairs(PairID,2); %Environment 1
        noise_std = 0.2; % Uniform Noise
        MemoryLimit=20;

        % Chosen neurons for visualization 
          Choose2NeuronsCA3=[1,2];% [6,10] and [8,53] . [6,10,15,18,19,38,66]

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

'Proportion Tagged Cells'
disp(cellfun(@(x) sum(x)/N, MemoryCells))

end

% Plot PreLoaded Data
PreLoadData=1;
MemoryStoreAntiCofiring_Diff;
