
lambda=0;

addpath("Data/")
load('/Path_To_/S3_C/Data/Results_Unsaturated.mat')

Activity=InputRates;
[DE_PV_In, neu_votesPV, D_Pos_PV] = plurality_voting_decoder(Activity(:,31:end), IntegerPos, 0, 1);

[DE_L, D_PosL] = linear_decoder(Activity(:,31:end), IntegerPos);
[DE_L1_In, sparsity,D_PosL1,TestPos] = linear_decoder_with_l1(Activity(:,31:end), IntegerPos, lambda);

figure;
scatter(IntegerPos(91:end), D_PosL1, 'b', 'filled');
hold on;
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');
title('Input Rates Regression');
xlabel('Expected Position');
ylabel('Predicted Position');

SaveDecoder{1,1}=D_PosL1';

scatter(IntegerPos(91:end), D_Pos_PV, 'k', 'filled');
hold on;
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');
title('Input Rates Regression');
xlabel('Expected Position');
ylabel('Predicted Position');
ylim([-20,50])
legend('Linear Decoder')


SaveDecoder{1,2}=D_Pos_PV;

'Linear Error Inputs'
DE_L1_In
'Assembly Tagging Inputs'
DE_PV_In



%%

num_iterations = 100; % Number of random subsamplings
num_neurons = size(Rates,1);
num_timepoints = size(IntegerPos,2) - 90; % Adjust for indexing

D_PosL_all = zeros(num_timepoints, num_iterations);
D_Pos_PV_all = zeros(num_timepoints, num_iterations);

num_example_trials = 5; % Number of example runs to plot
example_trials = randi(num_iterations, num_example_trials, 1); % Select random runs to overlay

for iter = 1:num_iterations
    % Randomly subsample neurons
    Activity = InputRates(randi(size(InputRates,1), num_neurons, 1)', :);
    
    % Decode using Linear Regression
    % [~, D_PosL] = linear_decoder(Activity(:,31:end), IntegerPos);
    [DE_L1, sparsity,D_PosL1,TestPos] = linear_decoder_with_l1(Activity(:,31:end), IntegerPos, lambda);

    % Decode using Plurality Voting
    [~, ~, D_Pos_PV] = plurality_voting_decoder(Activity(:,31:end), IntegerPos, 0, 1);
    
    % Store results
    D_PosL_all(:, iter) = D_PosL1;
    D_Pos_PV_all(:, iter) = D_Pos_PV;
end

% Compute mean and standard deviation
mean_D_PosL = mean(D_PosL_all, 2);
std_D_PosL = std(D_PosL_all, 0, 2);
mean_D_Pos_PV = mean(D_Pos_PV_all, 2);
std_D_Pos_PV = std(D_Pos_PV_all, 0, 2);

% Plot Regression Results with Shaded Error and Example Trials
figure; hold on;
scatter(IntegerPos(91:end), mean_D_PosL, 'b', 'filled');


SaveDecoder{2,1}=mean_D_PosL';


% fill([IntegerPos(91:end), fliplr(IntegerPos(91:end))], ...
%      [mean_D_PosL - std_D_PosL; flipud(mean_D_PosL + std_D_PosL)], ...
%      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');

% Overlay example trials
for i = 1:num_example_trials
    plot(IntegerPos(91:end), D_PosL_all(:, example_trials(i)), 'b', 'LineWidth', 0.8);
end

title('Regression');
xlabel('Expected Position');
ylabel('Predicted Position');

% Plot Plurality Voting Results with Shaded Error and Example Trials
hold on;
scatter(IntegerPos(91:end), mean_D_Pos_PV, 'k', 'filled');

SaveDecoder{2,2}=mean_D_Pos_PV';

% fill([IntegerPos(91:end), fliplr(IntegerPos(91:end))], ...
%      [mean_D_Pos_PV - std_D_Pos_PV; flipud(mean_D_Pos_PV + std_D_Pos_PV)], ...
%      'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');

% Overlay example trials
for i = 1:num_example_trials
    plot(IntegerPos(91:end), D_Pos_PV_all(:, example_trials(i)), 'k', 'LineWidth', 0.8);
end

title('Effect of Subsampling Inputs (5% Total Inputs = Total Outputs)');
xlabel('Expected Position');
ylabel('Predicted Position');
ylim([-20,50]);
legend('Linear Decoder')

load('/Path_To_/S3_C/Data/Results_Saturated.mat')

Activity=InputRates;
[DE_PV_In, neu_votesPV, D_Pos_PV] = plurality_voting_decoder(Activity(:,31:end), IntegerPos, 0, 1);

[DE_L, D_PosL] = linear_decoder(Activity(:,31:end), IntegerPos);
[DE_L1_In, sparsity,D_PosL1,TestPos] = linear_decoder_with_l1(Activity(:,31:end), IntegerPos, lambda);

figure;
scatter(IntegerPos(91:end), D_PosL1, 'b', 'filled');

SaveDecoder{3,1}=D_PosL1';

hold on;
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');
title('Saturated Inputs Regression');
xlabel('Expected Position');
ylabel('Predicted Position');


scatter(IntegerPos(91:end), D_Pos_PV, 'k', 'filled');

SaveDecoder{3,2}=D_Pos_PV;


hold on;
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');
title('Saturated Inputs Regression');
xlabel('Expected Position');
ylabel('Predicted Position');
ylim([-20,50])
legend('Linear Decoder')

'Linear Error Inputs'
DE_L1_In
'Assembly Tagging Inputs'
DE_PV_In


Activity=Rates;
[DE_PV_Out, neu_votesPV, D_Pos_PV] = plurality_voting_decoder(Activity(:,31:end), IntegerPos, 0, 1);

% lambda=0.05;

[DE_L1_Out, sparsity,D_PosL1,TestPos] = linear_decoder_with_l1(Activity(:,31:end), IntegerPos, lambda);
[DE_L, D_PosL] = linear_decoder(Activity(:,31:end), IntegerPos);

figure;
scatter(IntegerPos(91:end), D_PosL1, 'b', 'filled');

SaveDecoder{4,1}=D_PosL1';

hold on;
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');
title('Output Rates Assembly Tagging vs Linear Decoder');
xlabel('Expected Position');
ylabel('Predicted Position');


scatter(IntegerPos(91:end), D_Pos_PV, 'k', 'filled');

SaveDecoder{4,2}=D_Pos_PV;

hold on;
plot(min(IntegerPos(91:end)):max(IntegerPos(91:end)), min(IntegerPos(91:end)):max(IntegerPos(91:end)), 'r--');
title('Saturated Outputs Regression');
xlabel('Expected Position');
ylabel('Predicted Position');
ylim([-20,50])
legend('Linear Decoder')
'Assembly Tagging Inputs'
DE_L1_Out
'Assembly Tagging Outputs'
DE_PV_Out



SaveDecoders = vertcat(SaveDecoder{:}).';

