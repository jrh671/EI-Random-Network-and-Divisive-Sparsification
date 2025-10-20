addpath('./Network_Functions');

GenerateTuning=1;
GenerateActivity=1;
Plasticity=0;
PlotMaps=1;
VisualizePlasticity=0;
VisualizePositions=1;
Average=1;

% Parameters
N = 100;         % Number of neurons in CA3
K = 300;         % Number of inputs
P = 7;          % Number of positions per dimension (total positions = P x P)
T = 1200;         % Number of time steps
StartTime = 10; %End-time of visualization
EndTime = 100; %End-time of visualization
GoalRadius=0;
seed=1; %Environment 10
noise_std = 0; % Uniform Noise
MemoryLimit=20;

% Chosen neurons for visualization 
  Choose2NeuronsCA3=[2,7];% [6,10] and [8,53]

% Define special location and radius
special_location = [7, 2];
special_radius = 0; 

%% Run First Without Activity-Dependent-Plasticity
SetParameters4D;

if GenerateTuning==1
    Seed=seed;
    RunGenerateTuning;
end

if GenerateActivity==1
    'Begin Activity'
    RunGenerateActivity;
    FiringRates=FiringRates1;
    AllCells_RateMap;
    RateMap1=RateMap;
    nonZeroIndices1=nonZeroIndices;
end

Sort1=sorted_indices_CA3;

Visualize2DMaps;
SheetNoPlast =cumulative_activity_map_CA3(:,:,1);
Traj1=trajectory;
FRAve1=cumulative_activity_map_CA3;

CalcOverdispersion;

VisualizePlasticity=1;
QuickRun2;

DispersionStats;

VisualizePlasticity=0;
Visualize2DMaps;

SheetPlasticity=cumulative_activity_map_CA3(:,:,1);