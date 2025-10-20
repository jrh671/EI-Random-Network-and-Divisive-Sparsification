% --- Input ---
TC_raw = NS(neuron_indices, :)';  % size: bins × neurons
n_bins = size(TC_raw, 1);
n_neurons = size(TC_raw, 2);
angles = linspace(0, 2*pi, n_bins + 1);
angles(end) = [];

% --- Smoothing kernel (circular) ---
sigma = 1.5;  % adjust as needed
win_size = 5;  % kernel width in bins
x = -floor(win_size/2):floor(win_size/2);
gauss_kernel = exp(-x.^2 / (2 * sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);  % normalize

% --- Apply circular smoothing ---
TC = zeros(size(TC_raw));  % bins × neurons
pad_size = floor(length(gauss_kernel)/2);

for i = 1:n_neurons
    x = TC_raw(:, i);  % current tuning curve (length = n_bins)
    
    % Circular padding
    padded = [x(end-pad_size+1:end); x; x(1:pad_size)];
    
    % Convolve
    smoothed = conv(padded, gauss_kernel, 'same');
    
    % Extract smoothed output of original length
    TC(:, i) = smoothed(pad_size+1 : pad_size+n_bins);
end


% --- Normalize and clean ---
TC(TC < 0) = 0;  % enforce non-negativity

% --- Compute Rayleigh vector length ---
z = TC' * exp(1i * angles(:));
R = abs(z) ./ sum(TC, 1)';  % Rayleigh per neuron

% --- Compute Skaggs information (in bits) ---
P_i = ones(1, n_bins) / n_bins;  % uniform occupancy

FR_mean = sum(TC .* P_i', 1)';  % mean firing rate per neuron
epsilon = 1e-9;
TC_safe = TC + epsilon;
FR_mean_safe = FR_mean + epsilon;

Info = sum((TC_safe ./ FR_mean_safe') .* log2(TC_safe ./ FR_mean_safe') .* P_i', 1)';

% --- Additional condition: max value in NS > 0.3
max_vals = max(NS(neuron_indices, :), [], 2);
max_thresh = 0.4;
max_mask = max_vals > max_thresh;

% --- Apply thresholds ---
rayleigh_thresh = 0.9;
skaggs_thresh = 2;

% --- Combine all conditions
selected_mask = (R > rayleigh_thresh) & (Info > skaggs_thresh) & max_mask;
selected_indices = neuron_indices(selected_mask);

length(selected_indices)
% --- Select neuron subset
Neurons = selected_indices([2 5 13 15 14]);

% --- Choose subplot layout ---
subplot_shape = '2x2';  % options: '2x2' or '4x1'

if strcmp(subplot_shape, '2x2')
    subplot_rows = 2;
    subplot_cols = 2;
elseif strcmp(subplot_shape, '4x1')
    subplot_rows = 4;
    subplot_cols = 1;
else
    error('Invalid subplot_shape. Use "2x2" or "4x1".');
end

figure;

%% --- Subplot 1: neuron_strength_matrix (all + highlights, with gray background)
subplot(subplot_rows, subplot_cols, 1);
hold on;
plot(neuron_strength_matrix(neuron_indices, :)', 'Color', [0.8 0.8 0.8]);  % gray background
yline(4.9, 'k', 'LineWidth', 5);
highlight_colors = cell(1, length(Neurons));
for i = 1:length(Neurons)
    idx = Neurons(i);
    h = plot(neuron_strength_matrix(idx, :)', 'LineWidth', 5);  % highlighted line
    highlight_colors{i} = get(h, 'Color');  % store color
end
title('Feedforward Without Normalization');
xlabel('Stimulus (Position)')
ylabel('Linear Input Combination')
ylim([min(neuron_strength_matrix,[],'all'),max(neuron_strength_matrix,[],'all')])

%% --- Subplot 2: Thresholded neuron_strength_matrix (highlights only)
subplot(subplot_rows, subplot_cols, 2);
hold on;
threshold = 4.9;
for i = 1:length(Neurons)
    idx = Neurons(i);
    data = neuron_strength_matrix(idx, :);
    data(data < threshold) = 0;
    plot(data, 'LineWidth', 5, 'Color', highlight_colors{i});
end
title('Threshold Without Normalization');
xlabel('Stimulus (Position)')
ylabel('Threshold Linear Input Combination')

%% --- Subplot 3: NormM with background + highlights
subplot(subplot_rows, subplot_cols, 3);
hold on;
plot(NormM', 'Color', [0.8 0.8 0.8]);  % gray background
for i = 1:length(Neurons)
    idx = Neurons(i);
    plot(NormM(idx, :), 'LineWidth', 5, 'Color', highlight_colors{i});
end
yline(0.95, 'k', 'LineWidth', 5);
title('Feedforward With Normalization');
xlabel('Stimulus (Position)')
ylabel('Input Combination (Normalized)')
ylim([min(NormM,[],'all'),max(NormM,[],'all')])

% %% --- Subplot 4: NS tuning curves
% subplot(subplot_rows, subplot_cols, 4);
% hold on;
% for i = 1:length(Neurons)
%     idx = Neurons(i);
%     plot(NS(idx, :), 'LineWidth', 5, 'Color', highlight_colors{i});
% end
% title('Threshold With Normalization');
% xlabel('Stimulus (Position)')
% ylabel('Threshold Input Combination (Normalized)')
% 


%% --- Subplot 4: Thresholded neuron_strength_matrix (highlights only)
subplot(subplot_rows, subplot_cols, 4);
hold on;
threshold = 0.95;

for i = 1:length(Neurons)
    idx = Neurons(i);
    data = NormM(idx, :);
    data(data < threshold) = 0;
    plot(data, 'LineWidth', 5, 'Color', highlight_colors{i});
end

% yline(threshold, 'k', 'LineWidth', 5);
title('Thresholded Feedforward Inputs');
xlabel('Stimulus (Position)');
ylabel('Input Combination (Thresholded)');
ylim([0,1])








figure;
hold on;
threshold = 0.95;

for i = 1:length(Neurons)
    idx = Neurons(i);
    data = NormM(idx, :);
    data(data < threshold) = 0;
    plot(data, 'LineWidth', 5, 'Color', highlight_colors{i});
end

% yline(threshold, 'k', 'LineWidth', 5);
title('Thresholded Feedforward Inputs');
xlabel('Stimulus (Position)');
ylabel('Input Combination (Thresholded)');

data = NormM;
data(data < threshold) = 0;

figure;imagesc(data(neuron_indices,:));colormap bone


data = neuron_strength_matrix;
data(data < 4.9) = 0;
figure;imagesc(data(neuron_indices,:));colormap bone