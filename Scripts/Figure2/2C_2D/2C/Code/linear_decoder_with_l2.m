function [mean_squared_circular_error, sparsity, decoded_positions, positions_test] = linear_decoder_with_l2(Spikes, integer_pos, lambda)

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
    B = ridge(positions_train, Spikes_train', lambda, 0);
    B0 = B(1);
    B = B(2:end);

    decoded_positions = Spikes_test' * B + B0;

    % Calculate circular distance error
    num_positions = max(integer_pos);

    circular_distance = @(x, y, num_positions) min(abs(x - y), num_positions - abs(x - y));
    
    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), positions_test(i), num_positions);
        decoding_error = decoding_error + error^2;
    end
    
    % Mean squared circular error
    mean_squared_circular_error = decoding_error / length(decoded_positions);

    % Calculate sparsity s
    numNonZeroCoeffs = sum(B ~= 0);
    totalCoeffs = length(B);
    sparsity = (totalCoeffs - numNonZeroCoeffs) / totalCoeffs * 100;
end