load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Supplementary/S05_06/Data/1.mat')

% Number of samples (200 Time Bins Per second, 60 seconds)
Temp=200*60;

%Initialize Results
action_saveAvoid=zeros(25,Temp);
action_saveContext=zeros(25,Temp);

%Initialize Temporal Windows
last_three_seconds_activity=zeros(n_excit,3/mean_dt);
decision_time=ceil(5/mean_dt);
decision_interval=ceil(1/mean_dt);

for tt=1:Temp

% Record the spiking activity for contextual data
spiking_activity(:,tt) = spike_mat_excitContext(:,600*200+tt); 

% Update the last 3 seconds activity
last_three_seconds_activity = [last_three_seconds_activity(:,2:end), spiking_activity(:,tt)];

if tt == decision_time 
actions = RRN_DownstreamNeurons(last_three_seconds_activity,Indexical2); 
decision_time = decision_time + decision_interval;
action_saveAvoid(:,tt)=actions;

end
end

%Initialize Temporal Windows
last_three_seconds_activity=zeros(n_excit,3/mean_dt);
decision_time=ceil(5/mean_dt);
decision_interval=ceil(1/mean_dt);

for tt=1:Temp

% Record the spiking activity for Avoid context
spiking_activity(:,tt) = spike_mat_excitAvoid(:,tt); 

% Update the last 3 seconds activity
last_three_seconds_activity = [last_three_seconds_activity(:,2:end), spiking_activity(:,tt)];

% Run the Downstream Tagging Mechanism on the tagged neurons
if tt == decision_time
actions = RRN_DownstreamNeurons(last_three_seconds_activity,Indexical2); 
decision_time = decision_time + decision_interval;
action_saveContext(:,tt)=actions;

end

end

SumA=sum(action_saveAvoid);
SumC=sum(action_saveContext);

figure;subplot(1,2,1);stem(SumA,'g');title('Freeze Activations (Freeze Context)');ylim([0 5]);subplot(1,2,2);stem(SumC,'r');ylim([0 5]);title('Freeze Activations (Goal Context)');

Sheet1b=[SumA(:),SumC(:)];

sum(SumA)
sum(SumC)

figure; imagesc(FiringRateExb(Indexical2,31:330),[0 1])

Sheet2a=FiringRateExb(Indexical2,31:330);

FiringRateExAvoid = SlidingAverage_s(spike_mat_excitAvoid(:,600*200:end));

figure; imagesc(FiringRateExAvoid(Indexical2,:),[0 1])

Sheet2b=FiringRateExAvoid(Indexical2,:);
