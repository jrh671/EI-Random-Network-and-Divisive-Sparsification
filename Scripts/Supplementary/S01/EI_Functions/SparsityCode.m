
%% Compute Sparsity

% Compute the mean and standard error of Sparsity across instances (columns)
mean_Sparsity = mean(Sparsity, 2); % Mean across columns
stderr_Sparsity = std(Sparsity, 0, 2) / sqrt(size(Sparsity, 2)); % Standard error

if Sus==0
% Plot Sparsity with error bars
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
if Sus == 1
figure;
errorbar(linspace(0, 0.1, Iteration), mean_Susceptibility, stderr_Susceptibility, 'LineWidth', 1.5);
title('Susceptibility To Inhibition');
xlabel('Inhibitory Strength');
ylabel('Change in sparsity');
grid on;
end
