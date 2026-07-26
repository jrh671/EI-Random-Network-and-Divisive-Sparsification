PreLoad=1;

if PreLoad==1

load(['/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/' ...
    'Supplementary/S04/Conductance_Model/PreRun_Data/Fast_Pos/ConductanceModel_Data.mat'], ...
    'spike_mat_excit')

load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Supplementary/S04/Conductance_Model/PreRun_Data/Fast_Pos/ConductanceModel_Data.mat', 'positions')
end


% %% Compute Results
% pf_cell1=pf_cell;
% Run_EffectiveTuning;
[NS,Ordering]=Compute_RateMap(spike_mat_excit(:,1:end),positions(1:end));
mean_dt=total_time_vec(1);
neuron_indices=Ordering;
FiringRate=SlidingAverage_s(spike_mat_excit,mean_dt);
figure; imagesc(FiringRate(neuron_indices,:),[0 1]);
title('Neuron Activity');xlabel('Time (s)'); ylabel('Neuron ID (Sorted By Position)')
n_pos=30;
% % Assembly Voting Decoder
IntegerPos=ceil(SlidingAverage_s_mean(positions,mean_dt));
time_bin_length=1; excluded_p=0;
[DE_PVo, neu_votesPV, D_Pos_PVo] = plurality_voting_decoder(FiringRate(:,1:end), IntegerPos, excluded_p, time_bin_length);

figure;
plot(IntegerPos(length(IntegerPos)-length(D_Pos_PVo)+1:end));hold on; plot(D_Pos_PVo, 'bo')
xlabel('Time Point');ylabel("Location"); title('Decoding Position From Activity');legend('True Location','Predicted Location')


PosAveSheet=NS(Ordering,:);
FRSheet=FiringRate(neuron_indices,:);
DecodeSheet_True=IntegerPos(length(IntegerPos)-length(D_Pos_PVo)+1:end)';
DecodeSheet_Pred=D_Pos_PVo';