addpath('SavedFiles_Fast/')
addpath('Helper_Functions/')

load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1G/PreRun_Data/PF_CellFastSTDP1.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1G/PreRun_Data/resultsFastSTDP1.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1G/PreRun_Data/W_InputEFastSTDP1.mat')


threshold=0.92;
n_pos=30;
n_laps=5;
Run_EffectiveTuning
figure;imagesc(results(neuron_indices,:), [ 0 1])
SaveRates=results(neuron_indices,:);