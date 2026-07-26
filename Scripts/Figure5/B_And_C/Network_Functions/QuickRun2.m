
%% Run Next With Activity-Dependent-Plasticity
Plasticity=1;

SetParameters4C;
if GenerateTuning==1
    RunGenerateTuning;
end

if GenerateActivity==1
    'Begin Activity'
    RunGenerateActivity;
end

FiringRates=FiringRates2;
Traj2=trajectory;
FRAve2=cumulative_activity_map_CA3;
AllCells_RateMap;
RateMap2=RateMap;

% Visualization
Visualize2DMaps;



