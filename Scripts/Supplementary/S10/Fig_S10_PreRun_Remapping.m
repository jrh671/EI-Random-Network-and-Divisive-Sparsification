addpath('./Network_Functions');
addpath('./PreRun_Data');

AAB=0; %1 = No intermediate remapping. 0 = With intermediate remapping
PairID=1; %Examples used in Figure SX: 1,8,11,19,52,55

if AAB==1
    load('./PreRun_Data/StoreMemory_AAB.mat', 'TopPairs_AAB'); 
    TopPairs=TopPairs_AAB;
else
    load('./PreRun_Data/StoreMemory_ABA.mat', 'TopPairs_ABA'); 
    TopPairs=TopPairs_ABA;
end

Single=1;

GenerateTuning=1;
GenerateActivity=1;
Plasticity=1;
PlotMaps=1;
VisualizePlasticity=0;
VisualizePositions=0;
PlotAntiCo=0;
PlotStability=0;
Rates = cell(8,1);
Sorting = cell(8,1);
MemoryCells = cell(8,1);
TemporalRates = cell(8,1);
TemporalSpikes = cell(8,1);


ContextGoal = 1; 



for Context = ContextGoal
    clear MCells MemoryCells

    if Context == 1
        Epochs = 1:8;
    elseif Context == 2
        Epochs = 1:8;
    end


    for Epoch = Epochs
    
        fprintf('Epoch %d\n', Epoch);
    
        % Parameters
        N = 100;         % Number of neurons in CA3
        K = 300;         % Number of inputs
        P = 10;          % Number of positions per dimension (total positions = P x P)
        T = 1200;         % Number of time steps
        StartTime = 50; %End-time of visualization
        EndTime = 150; %End-time of visualization
        GoalRadius=0;
        seed1=TopPairs(PairID,1); %Environment 0
        seed2=TopPairs(PairID,2); %Environment 1
        noise_std = 0.2; % Uniform Noise
        MemoryLimit=20;
        % Chosen neurons for visualization 
          Choose2NeuronsCA3=[6,7];
        
        % Define special location and radius
        special_location = [7, 2];
    


        if Epoch == 1 
            special_radius = GoalRadius; 
        else
            special_radius = 0;
        end

       

    
        
        %% Run First Without Activity-Dependent-Plasticity
        
        SetParameters4D;



        if Epoch==1
            if GenerateTuning==1
                Plasticity=1;
                Seed = seed1;
                RunGenerateTuning;
                SaveConnections1 = A_CA3;
                SaveTuning1 = B;
            end

        elseif Epoch == 2
            
            if GenerateTuning==1
                Plasticity=0;
                Seed = seed1;
                RunGenerateTuning;
                SaveConnections1 = A_CA3;
                SaveTuning1 = B;
            end

        elseif Epoch == 3 
            Plasticity=0;
                Seed = seed2;
                Remap_InputTuning;
                SaveConnections2 = A_CA3;
                SaveTuning2 = B;
        elseif Epoch == 4
                Plasticity=0;
                Seed = seed1;
                Remap_InputTuning;
                SaveTuning3 = B;
        elseif Epoch == 5
                Plasticity=0;
                Seed = seed2;
                Remap_InputTuning;
                SaveTuning3 = B;
        elseif Epoch == 6
                Plasticity=1;
                Seed = seed1;
                Remap_InputTuning;
                SaveTuning3 = B;
        elseif Epoch == 7
                Plasticity=1;
                Seed = seed2;
                Remap_InputTuning;
                SaveTuning3 = B;
        elseif Epoch == 8
                Plasticity=1;
                Seed = seed2;
                Remap_InputTuning;
                SaveTuning3 = B;
        end
        
        if GenerateActivity==1
            RunGenerateActivity;
            if Epoch>0
            Rates{Epoch} = firing_rate_snapshots_CA3;
            Sorting{Epoch} = sorted_indices_CA3;
            MemoryCells{Epoch} = memory_neurons_CA3;
            end
        end
    
    end

    if Epoch>0
    Epochs = 1:8;
    Visualize2DMaps;
    end

'Proportion Tagged Cells'
disp(cellfun(@(x) sum(x)/N, MemoryCells))

end

