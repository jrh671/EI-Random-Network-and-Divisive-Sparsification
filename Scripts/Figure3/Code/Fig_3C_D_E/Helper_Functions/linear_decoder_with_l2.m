function [mean_squared_circular_error, sparsity, decoded_positions, positions_test] = linear_decoder_with_l2(Spikes, integer_pos, lambda)
    % Function to calculate circular distance
    circular_distance = @(x, y, num_positions) min(abs(x - y), num_positions - abs(x - y));

    % Split data into training and testing sets
    num_samples = size(Spikes, 2);
    half_samples = floor(num_samples / 2);

    % Training data (first half)
    Spikes_train = Spikes(:, 1:half_samples);
    positions_train = integer_pos(1:half_samples)';  % Transpose to make it a column vector [samples, 1]

    % Testing data (second half)
    Spikes_test = Spikes(:, half_samples+1:end);
    positions_test = integer_pos(half_samples+1:end);

    % Train Ridge regression model
    % Ridge requires predictors (Spikes_train') with samples as rows and features as columns
    % Lambda is the regularization parameter (L2 regularization strength)
    B = ridge(positions_train, Spikes_train', lambda, 0);

    % Ensure B has the same number of elements as the features
    % Trim the extra intercept term if present
    if size(B, 1) > size(Spikes_train, 1)
        B = B(2:end);  % Remove the first element, which is likely the bias term
    end

    % In ridge regression, the intercept must be calculated manually
    % Compute intercept (B0) using the means of training data
    mean_train_features = mean(Spikes_train, 2);  % Column vector of [num_features, 1]
    B0 = mean(positions_train) - sum(mean_train_features .* B);

    % Decode positions using the Ridge model on the test set
    decoded_positions = Spikes_test' * B + B0;

    % Calculate circular distance error
    num_positions = max(integer_pos);

    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), positions_test(i), num_positions);
        decoding_error = decoding_error + error^2;
    end

    % Mean squared circular error
    mean_squared_circular_error = decoding_error / length(decoded_positions);

    % Calculate sparsity (Ridge generally does not force coefficients to zero, so sparsity is usually 0)
    numNonZeroCoeffs = sum(B ~= 0);
    totalCoeffs = length(B);
    sparsity = (totalCoeffs - numNonZeroCoeffs) / totalCoeffs * 100;
end
