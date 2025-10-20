%% Compute Decoding Error

% Calculating the mean and standard error for each condition
means = mean(data, 2); % Mean across the second dimension
std_errors = std(data, [], 2) / sqrt(size(data, 2)); % Standard Error of the Mean

% Conditions (for x-axis)
if Across==1
conditions = linspace(0,0.1,Iteration); % Adjust this as per your conditions
Colors = 'k';
elseif Across==2
conditions1 = 0.0500*ones(1,Iteration);
conditions2 = linspace(0,0.4,Iteration); % Adjust this as per your conditions
elseif Across==3
conditions = linspace(0,0.4,Iteration); % Adjust this as per your conditions
Colors = 'b';
end

if Instance ==1
figure

end

Colors=['g','m','c','y','r'];
hold on;
hErrorbar = errorbar(conditions(1:Iteration), means(1:Iteration), std_errors(1:Iteration), '-', 'Color',  Colors(Instance), 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
% ylim([0,1.1])
title('Predicted vs Actual Position MSE');