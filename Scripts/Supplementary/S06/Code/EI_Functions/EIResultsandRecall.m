clear FiringRateExa
    
FiringRateExa = SlidingAverage_s(spike_mat_excitf);
CollectAvoidRates{J}=FiringRateExa;

%% Place Avoidance Training and Tagging

TrainGoal;
TaggedSize(J)=NumTagged;

pf_cell=pf_cell1;
Run_EffectiveTuning;
if J==1
figure;imagesc(FiringRateExa(neuron_indices,:),[0 1]);
end
%% Contextual Freezing Training and Tagging

% pf_cell=pf_cell1;
%     TrainingPrimacyCells;
% spike_mat_excitAvoid=spike_mat_excitB;

pf_cell=pf_cell2;
    TrainingPrimacyCells;
spike_mat_excitContext=spike_mat_excitB;

%% Re-expose to both context and look for reactivation (Freeze and Avoidance responses)
if J==1
figure;imagesc(FiringRateExb(Indexical2,:),[0 1]);
end

if Recall == 1
[SpikesRecall,positionsRecall,~,~,ReverseData,FreezeData, ~,~,FR] = RRN_SensorimotorRecallNetwork(n_input,n_excit,n_inhib,W_inputE,pf_cell1,pf_cell2,TopNeurons,Indexical2);
 
CollectRates{J}=FR;
PosRecall=SlidingAverage_s_mean(positionsRecall,.005);

pf_cell=pf_cell1;
    Run_EffectiveTuning;

CollectIndices{1,J} = neuron_indices;

pf_cell=pf_cell2;
    Run_EffectiveTuning;

CollectIndices{2,J} = neuron_indices;

if J==1

figure;imagesc(FR(CollectIndices{1,1},:),[0 1]);colorbar


figure;imagesc(FR(CollectIndices{2,1},:),[0 1]);colorbar
end
end