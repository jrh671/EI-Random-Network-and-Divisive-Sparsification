function [mean_squared_circular_error, sparsity, decoded_positions, positions_test] = linear_decoder_with_l1(Spikes, integer_pos, lambda)
    % Function to calculate circular distance
    circular_distance = @(x, y, num_positions) min(abs(x - y), num_positions - abs(x - y));

    % Split data into training and testing sets
    num_samples = size(Spikes, 2);
    half_samples = floor(num_samples / 2);

    % Training data (first half)
    Spikes_train = Spikes(:, 1:half_samples);
    positions_train = integer_pos(1:half_samples);

    % Testing data (second half)
    Spikes_test = Spikes(:, half_samples+1:end);
    positions_test = integer_pos(half_samples+1:end);

    % Train Lasso regression model
    % [B, FitInfo] = lasso(Spikes_train', positions_train, 'Lambda', lambda);
    [B, FitInfo] = lasso(Spikes_train', positions_train, 'Lambda', lambda);

    % Select the best model
    B0 = FitInfo.Intercept;
    Coefficients = B(:, 1);

    % Decode positions using the lasso model on the test set
    decoded_positions = round(Spikes_test' * Coefficients + B0);

    % Calculate circular distance error
    num_positions = max(integer_pos);

    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), positions_test(i), num_positions);
        decoding_error = decoding_error + error^2;
    end

    % Mean squared circular error
    mean_squared_circular_error = decoding_error / length(decoded_positions);

    % Calculate sparsity
    numNonZeroCoeffs = sum(Coefficients ~= 0);
    totalCoeffs = length(Coefficients);
    sparsity = (totalCoeffs - numNonZeroCoeffs) / totalCoeffs * 100;
end
