addpath('./Parallel');
addpath('./PreRun_Data');
addpath('./Helper_Functions');

n_pos=30;
n_laps=5;
integer_pos=mod(1:n_pos*n_laps,n_pos);
mean_dt=0.0005; 

load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/resultsFast.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/W_InputEFast.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/PF_CellFast.mat')

TM=1;
EI_FullAnalysis_Helper;

load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/resultsSlow.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/W_InputESlow.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure1/Code/1H_1I/PreRun_Data/PF_CellSlow.mat')

TM=2;
mean_dt=0.005; 
EI_FullAnalysis_Helper;

EI_Rates_PreRun;