function [normalized_matrix, sorted_matrix, sort_order] = compute_effective_tuning(input_neurons, weight_matrix, n_pos, threshold)
    % Input:
    % input_neurons: 1x1000 cell array, each cell containing a vector of positions (1 to 30)
    % weight_matrix: 1000x500 matrix of connection strengths between input neurons and recurrent network neurons
    % n_pos: Number of positions (30)
    % threshold: Minimum normalized value (0 to 1) to include rows in the sorting
    
    % Output:
    % normalized_matrix: Column-wise normalized neuron strength matrix
    % sorted_matrix: Row-sorted matrix based on the column of the highest value
    % sort_order: Row indices corresponding to the sorting order
    
    % Step 1: Convert input_neurons into a binary matrix
    input_matrix = zeros(1000, n_pos);
    for i = 1:length(input_neurons)
        input_matrix(i, ceil(input_neurons{i})) = 1;
    end
    
    % Step 2: Compute the neuron strength matrix
    neuron_strength_matrix = (input_matrix' * weight_matrix)'; 
    
    % Step 3: Normalize neuron_strength_matrix column-wise by maximum row value in each column
    max_vals = max(neuron_strength_matrix, [], 1);
    normalized_matrix = neuron_strength_matrix ./ max_vals;
    
    % Step 4: Determine rows above the threshold
    max_row_vals = max(normalized_matrix, [], 2); % Find the max value for each row
    above_threshold = max_row_vals >= threshold; % Boolean array for rows meeting the threshold
    below_threshold = ~above_threshold; % Boolean array for rows below the threshold
    
    % Step 5: Sort rows above the threshold by their highest value's column
    [~, max_col_indices] = max(normalized_matrix(above_threshold, :), [], 2);
    [~, sort_order_above] = sort(max_col_indices);
    above_sorted_rows = find(above_threshold); % Indices of rows above the threshold
    sorted_above = above_sorted_rows(sort_order_above); % Sorted rows above the threshold
    
    % Step 6: Append rows below the threshold at the end in arbitrary order
    below_sorted_rows = find(below_threshold); % Indices of rows below the threshold
    sort_order = [sorted_above; below_sorted_rows]; % Final sorting order
    
    % Step 7: Apply the sorting order to the matrix
    sorted_matrix = normalized_matrix(sort_order, :);
end
