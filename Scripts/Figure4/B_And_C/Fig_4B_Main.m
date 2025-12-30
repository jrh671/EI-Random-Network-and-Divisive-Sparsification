%% EI Random NETWORK (Main)
addpath('./Network_Functions');
addpath('./PreRun_Data');
load('./PreRun_Data/TwoPFs.mat')

GenerateTuning=1;
GenerateActivity=1;
VisualizePlasticity=0;
Plasticity=0;

% Parameters
N = 100;         % Number of neurons in CA3
K = 300;         % Number of inputs
P = 7;          % Number of positions per dimension (total positions = P x P)
T = 5000;         % Number of time steps
noise_std=0;
MemoryLimit=1000;

% Chosen neurons for visualization 
Choose2NeuronsCA3=[59,68]; % 2D Compatible Neurons: [7,16,93,99,97,100]]

% Define special location and radius
special_location = [7, 2];
special_radius = 0; 

SetParameters4B;

seed1=0;% Random Seeds (Map 1)
seed2=1; % Random Seeds (Map 2)

% Generative Model
if GenerateTuning==1
    Seed=seed1;
    RunGenerateTuning;
end

% Generate Activity
if GenerateActivity==1
    RunGenerateActivity;
    SaveConnections1 = A_CA3;
    SaveTuning1 = B;
end

% Visualization
Visualize2DMaps;
% 
PlaceCellStats;

SmallPFs1=find(numPeaks>0);
commonValues1 = intersect(TwoPFs, SmallPFs1);

% Remapping Function
Seed=seed2;
Remap_InputTuning;
SaveConnections2 = A_CA3;
SaveTuning2 = B;

% Generate Activity
RunGenerateActivity;

% Visualization
Visualize2DMaps;

PlaceCellStats;

SmallPFs2=find(numPeaks>0);
commonValues2 = intersect(TwoPFs, SmallPFs2);
