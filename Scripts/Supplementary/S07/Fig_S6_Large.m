addpath('./Network_Functions');

GenerateTuning=1;
GenerateActivity=1;
Plasticity=0;
VisualizePlasticity=0;

Small=0; %Do Not Change (For Recording Data)
Large=1; %Do Not Change (For Recording Data)

% Parameters
N = 100;         % Number of neurons in CA3
K = 300;         % Number of inputs
P = 20;          % Number of positions per dimension (total positions = P x P)
T = 8000;         % Number of time steps
EndTime = 300; % End-time of visualization
seed=1;  % Environment ID. 
noise_std=0;
MemoryLimit=1000; 

% For Other Seeds choose different Neurons Below:
Choose2NeuronsCA3= [20,34];%[14,62];%

% Define special location and radius
special_location = [7, 2]; 
special_radius = 0; 

%% Run First Without Activity-Dependent-Plasticity
SetParametersLarge;

if GenerateTuning==1
    Seed=seed;
    RunGenerateTuning;
end

if GenerateActivity==1
    'Begin Activity'
    RunGenerateActivity;
end


Visualize2DMaps;

PlaceCellStats;

seed=2;  % Environment ID. 

%% Run Next With Activity-Dependent-Plasticity
Plasticity=1;

if GenerateTuning==1
    Seed=seed;
    Remap_InputTuning;
end

if GenerateActivity==1
    'Begin Activity'
    RunGenerateActivity;
end

% Visualization
Visualize2DMaps;

LargePFs=find(numPeaks>1 & numPeaks<4);


