%% Compute Results
pf_cell1=pf_cell;
if Iteration == 1
Run_EffectiveTuning;
end
FiringRate=SlidingAverage_s(spike_mat_excit,mean_dt);
if Instance == 1 && Across == 1 && Option == 1
figure; imagesc(FiringRate(neuron_indices,n_pos+1:end),[0 1]);title('Firing Rates');xlabel('Time');ylabel('Neurons')
end