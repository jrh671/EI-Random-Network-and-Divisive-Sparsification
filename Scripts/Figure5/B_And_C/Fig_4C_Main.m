addpath('./Network_Functions');

GenerateTuning=1;
GenerateActivity=1;
Plasticity=0;
VisualizePlasticity=0;

% Parameters
N = 100;         % Number of neurons
K = 300;         % Number of inputs
P = 20;          % Number of positions per dimension (total positions = P x P)
T = 8000;         % Number of time steps
EndTime = 300; % End-time of visualization
seed=1;  % Environment ID. 
noise_std=0;
MemoryLimit=1000; 

% For Visualization
Choose2NeuronsCA3= [59,68];

% Define special location and radius
special_location = [7, 2]; 
special_radius = 0; 

%% Run First Without Activity-Dependent-Plasticity
SetParameters4C;

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
    Traj1=Trajectory1;
    nonZeroIndices1=nonZeroIndices;
end


Visualize2DMaps;
PlaceCellStats;



%% Run Next With Activity-Dependent-Plasticity
Plasticity=1;

SetParameters4C;
if GenerateTuning==1
    RunGenerateTuning;
end

if GenerateActivity==1
    'Begin Activity'
    RunGenerateActivity;
    FiringRates=FiringRates2;
    AllCells_RateMap;
    RateMap2=RateMap;
    Traj2=Trajectory2;
    nonZeroIndices2=nonZeroIndices;

end

% Visualization
Visualize2DMaps;
TwoPFs=find(numPeaks==2);
Overdispersion