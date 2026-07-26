%% Run Compute Relative Strength
input_neurons = pf_cell; % 1x1000 cell array, each cell contains a vector of positions
weight_matrix = W_inputE; % 1000x500 matrix with random connection strengths

% Compute the relative strengths and neuron indices
[NormM, SortM, neuron_indices] = compute_effective_tuning(input_neurons, weight_matrix, n_pos, threshold);
