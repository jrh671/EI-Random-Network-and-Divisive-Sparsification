addpath('SavedFiles_Fast/')
addpath('Helper_Functions/')

load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1H/PreRun_Data/PF_CellFastBCM1.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1H/PreRun_Data/resultsFastBCM1.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1H/PreRun_Data/W_InputEFastBCM1.mat')


threshold=0;
n_pos=30;
n_laps=5;
Run_EffectiveTuning
figure;imagesc(results(neuron_indices,:), [ 0 1])
SaveRates=results(neuron_indices,:);