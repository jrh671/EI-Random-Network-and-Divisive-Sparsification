addpath('./Parallel');
addpath('./PreRun_Data');


for TM=1:2

if TM==1

load('PF_CellSlow.mat')
load('resultsSlow.mat')
load('W_InputESlow.mat')
else

load('PF_CellFast.mat')
load('resultsFast.mat')
load('W_InputEFast.mat')

end

n_pos=30;
n_laps=5;
integer_pos=mod(1:n_pos*n_laps,n_pos);


for idx = 1:length(results)

Run_Decoder3A; %In Helper_Functions Folder

end

if TM==1
    Data1=data;
    ChooseData=Data1;
else
    Data2=data;
    ChooseData=Data2;
end


%% Compute Decoding Error
data=ChooseData;

% Calculating the mean and standard error for each condition
means = mean(data, 2); % Mean across the second dimension
std_errors = std(data, [], 2) / sqrt(size(data, 2)); % Standard Error of the Mean

% Conditions (for x-axis)
conditions = linspace(0,0.3,Iteration); % Adjust this as per your conditions
Colors = ['r', 'b'];


hErrorbar = errorbar(conditions(1:Iteration), means(1:Iteration), std_errors(1:Iteration), 'o', 'Color',  Colors(TM), 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
hold on;
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
ylim([0,85])
title('Predicted vs Actual Position MSE');

SaveMeans{TM}=means(1:Iteration);
SaveSEM{TM}=std_errors(1:Iteration);

end