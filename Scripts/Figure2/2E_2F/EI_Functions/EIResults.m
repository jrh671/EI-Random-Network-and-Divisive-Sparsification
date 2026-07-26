%% Compute Results
pf_cell1=pf_cell;
Run_EffectiveTuning;
[NS,Ordering]=Compute_RateMap(spike_mat_excit(:,(n_pos/mean_dt)+1:end),positions((n_pos/mean_dt)+1:end));

FiringRate=SlidingAverage_s(spike_mat_excit,mean_dt);
figure; imagesc(FiringRate(neuron_indices,n_pos+1:end),[0 1]);
title('Neuron Activity');xlabel('Time (s)'); ylabel('Neuron ID (Sorted By Position)')

% Assembly Voting Decoder
IntegerPos=SlidingAverage_s_mean(integer_pos,mean_dt);
time_bin_length=1;excluded_p=0;
[DE_PVo, neu_votesPV, D_Pos_PVo] = plurality_voting_decoder(FiringRate(:,n_pos+1:end), IntegerPos, excluded_p, time_bin_length);

figure;
plot(IntegerPos(length(IntegerPos)-length(D_Pos_PVo)+1:end));hold on; plot(D_Pos_PVo, 'bo')
xlabel('Time Point');ylabel("Location"); title('Decoding Position From Activity');legend('True Location','Predicted Location')

figure;scatter(IntegerPos(length(IntegerPos)-length(D_Pos_PVo)+1:end),D_Pos_PVo);hold on;
hold on;plot(linspace(0,n_pos,n_pos),linspace(0,n_pos,n_pos));xlabel('Predicted Location');ylabel("Actual Location"); title('Decoding Position From Activity')