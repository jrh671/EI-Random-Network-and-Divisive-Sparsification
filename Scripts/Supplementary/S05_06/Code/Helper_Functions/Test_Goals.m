load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S05_06/Data/1.mat')

TrainGoal;

% Number of samples (200 Time Bins Per second, 60 seconds)
Temp=200*60;

%Initialize Results
action_saveAvoid=zeros(1,Temp);
action_saveContext=zeros(1,Temp);

%Initialize Temporal Windows
last_three_seconds_activity=zeros(n_excit,3/mean_dt);
decision_time=ceil(5/mean_dt);
decision_interval=ceil(1/mean_dt);

for tt=1:Temp
    
% Record the spiking activity for contextual data
spiking_activity(:,tt) = spike_mat_excitContext(:,tt); % Update this with your own code

% Update the last 3 seconds activity
last_three_seconds_activity = [last_three_seconds_activity(:,2:end), spiking_activity(:,tt)];

% Run the Downstream Tagging Mechanism on the tagged neurons
if tt == decision_time 
actions = RRN_DownstreamNeuron(last_three_seconds_activity,TopNeurons); 
decision_time = decision_time + decision_interval;
action_saveContext(:,tt)=actions;
end

end


last_three_seconds_activity=zeros(n_excit,3/mean_dt);
decision_time=ceil(5/mean_dt);
decision_interval=ceil(1/mean_dt);

for tt=1:Temp
% Record the spiking activity for avoidance data
spiking_activity(:,tt) = spike_mat_excitAvoid(:,tt); 

% Update the last 3 seconds activity
last_three_seconds_activity = [last_three_seconds_activity(:,2:end), spiking_activity(:,tt)];

% Run the Downstream Tagging Mechanism on the tagging-associated neuron
if tt == decision_time
actions = RRN_DownstreamNeuron(last_three_seconds_activity,TopNeurons); 
decision_time = decision_time + decision_interval;
action_saveAvoid(:,tt)=actions;
end
end

%% Plot Results
SumA=action_saveAvoid;
SumC=action_saveContext;

figure;subplot(1,2,1);stem(SumA,'g');title('Goal Activations (Goal Context)');ylim([0,2]);subplot(1,2,2);stem(SumC,'r');title('Goal Activations (Freeze Context)');ylim([0,2]);

Sheet3b=[SumA(:),SumC(:)];

sum(SumA)
sum(SumC)