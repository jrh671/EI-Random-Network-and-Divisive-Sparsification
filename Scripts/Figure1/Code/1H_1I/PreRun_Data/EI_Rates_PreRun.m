    
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/resultsFast.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/W_InputEFast.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/PF_CellFast.mat')

pf_cell=pf_cellC;
W_inputE=W_inputECell{1};

n_pos=30;n_laps=5;
integer_pos=mod(1:n_pos*n_laps,n_pos);
threshold=0.99;
Run_EffectiveTuning;

Option=1;
Iterations=[1,7,14];

Data1=results;


for I=1:3
    
Iteration=Iterations(I);

figure
    % Compute idx
idx = (Iteration - 1) * 10 + Option;
imagesc(results{idx}(neuron_indices,:),[0 1])
SaveFR_Fast{I}=results{idx}(neuron_indices,:);

title((Iteration-1)*0.015)
xlabel('Time')
ylabel('Neuron (Sorted By Position)')
pause(1)
end

load('/Users/josehurtado/Documents/MATLAB/Manuscript/Home_Dir/Random_EIPlaceNet/Figure2/Code/Fig_2A_B/PreRun_Data/resultsSlow.mat')
load('/Users/josehurtado/Documents/MATLAB/Manuscript/Home_Dir/Random_EIPlaceNet/Figure2/Code/Fig_2A_B/PreRun_Data/W_InputESlow.mat')
load('/Users/josehurtado/Documents/MATLAB/Manuscript/Home_Dir/Random_EIPlaceNet/Figure2/Code/Fig_2A_B/PreRun_Data/PF_CellSlow.mat')

pf_cell=pf_cellC;
W_inputE=W_inputECell{1};

pf_cell=pf_cellC;
threshold=0.97;
Run_EffectiveTuning;

Data2=results;

Option=1;
for I=1:3
figure

Iteration=Iterations(I);
% Compute idx
idx = (Iteration - 1) * 10 + Option;
imagesc(results{idx}(neuron_indices,:),[0 1])
title((Iteration-1)*0.015)
xlabel('Time')
ylabel('Neuron (Sorted By Position)')
pause(1)

SaveFR_Slow{I}=results{idx}(neuron_indices,:);
end