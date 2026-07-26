load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1C_1D/PreRun_Data/Data.mat')

Iteration=21;
%% Compute Sparsity
% Compute the mean and standard error of Sparsity across instances (columns)
mean_Sparsity = mean(Sparsity, 2); % Mean across columns
stderr_Sparsity = std(Sparsity, 0, 2) / sqrt(size(Sparsity, 2)); % Standard error

% Plot Sparsity with error bars
Instance=1;
if Instance == 1
figure;
errorbar(linspace(0, 0.1, Iteration), mean_Sparsity, stderr_Sparsity, 'LineWidth', 1.5);
title('Sparsity (Proportion Inactive Neurons)');
xlabel('Inhibition');
ylabel('Sparsity');
grid on;
end
% Compute the gradient (difference of adjacent elements) for each instance
Susceptibility = diff(Sparsity);

% Append a zero at the start for each column
Susceptibility = [zeros(1, size(Susceptibility, 2)); Susceptibility];

% Compute the mean and standard error of Susceptibility across instances
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



load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure1/Code/1C_1D/PreRun_Data/Single_FR.mat')

for i = 5:4:21
figure;
imagesc(FR{i},[0 1]);title('Firing Rates');xlabel('Time');ylabel('Neurons')
end