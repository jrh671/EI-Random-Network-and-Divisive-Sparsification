load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1C_1D/PreRun_Data/Data.mat')

Iteration=21;
%% Compute Sparsity
mean_Sparsity = mean(Sparsity, 2);
stderr_Sparsity = std(Sparsity, 0, 2) / sqrt(size(Sparsity, 2)); % Standard error

Instance=1;
if Instance == 1
figure;
errorbar(linspace(0, 0.1, Iteration), mean_Sparsity, stderr_Sparsity, 'LineWidth', 1.5);
title('Sparsity (Proportion Inactive Neurons)');
xlabel('Inhibition');
ylabel('Sparsity');
grid on;
end

Susceptibility = diff(Sparsity);

Susceptibility = [zeros(1, size(Susceptibility, 2)); Susceptibility];

mean_Susceptibility = mean(Susceptibility, 2);
stderr_Susceptibility = std(Susceptibility, 0, 2) / sqrt(size(Susceptibility, 2));

% Plot Susceptibility with error bars
if Instance == 1
figure;
errorbar(linspace(0, 0.1, Iteration), mean_Susceptibility, stderr_Susceptibility, 'LineWidth', 1.5);
title('Susceptibility To Inhibition');
xlabel('Inhibitory Strength');
ylabel('Change in sparsity');
grid on;
end

%% Compute Decoding Error

for Across = 1:3
    Instance=Across;
data=Data{Across};

means = mean(data, 2); 
std_errors = std(data, [], 2) / sqrt(size(data, 2)); 

SaveMeans{Across}=means;
SaveSTDs{Across}=std_errors;


% Conditions (for x-axis)
if Across==1
conditions = linspace(0,0.1,Iteration); 
Colors = 'k';
elseif Across==2
conditions1 = 0.0500*ones(1,Iteration);
conditions2 = linspace(0,0.4,Iteration); 
elseif Across==3
conditions = linspace(0,0.4,Iteration); 
Colors = 'b';
end

if Instance ==1
figure
hold on;
hErrorbar = errorbar(conditions(1:Iteration), means(1:Iteration)/max(means(1:Iteration)), std_errors(1:Iteration)/max(means(1:Iteration)), 'o', 'Color',  Colors, 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
ylim([0,1.1])
title('Predicted vs Actual Position MSE');
end

if Instance == 2
Colors = 'r';
hold on;
hErrorbar = errorbar(conditions1(1:Iteration), means(1:Iteration)/max(means(1:Iteration)), std_errors(1:Iteration)/max(means(1:Iteration)), 'o', 'Color',  Colors, 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
ylim([0,1.1])
title('Predicted vs Actual Position MSE');
    
figure
Colors = 'r';
hold on;
hErrorbar = errorbar(conditions2(1:Iteration), means(1:Iteration)/max(means(1:Iteration)), std_errors(1:Iteration)/max(means(1:Iteration)), 'o', 'Color',  Colors, 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
ylim([0,1.1])
title('Predicted vs Actual Position MSE');
end

if Instance == 3
hold on;
hErrorbar = errorbar(conditions(1:Iteration), means(1:Iteration)/max(means(1:Iteration)), std_errors(1:Iteration)/max(means(1:Iteration)), 'o', 'Color',  Colors, 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
ylim([0,1.1])
title('Predicted vs Actual Position MSE');
end

end

