addpath('SavedFiles_Fast/')
addpath('Helper_Functions/')

load('/Directory_To/SavedFiles_Fast/PF_CellFastSTDP1.mat')
load('/Directory_To/SavedFiles_Fast/resultsFastSTDP1.mat')
load('/Directory_To/SavedFiles_Fast/W_InputEFastSTDP1.mat')


threshold=0.975;
n_pos=30;
n_laps=5;
Run_EffectiveTuning
figure;imagesc(results(neuron_indices,:), [ 0 1])
SaveRates=results(neuron_indices,:);