function mean_squared_circular_error = linear_decoder2(Spikes, integer_pos)
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

    % Center the training data and labels
    Spikes_mean = mean(Spikes_train, 2);  % Mean across columns (features)
    positions_mean = mean(positions_train);  % Mean of training labels

    % Subtract mean to center the training data and labels
    Spikes_train_centered = Spikes_train - Spikes_mean;
    positions_train_centered = positions_train - positions_mean;

    % Center the testing data using training mean (important!)
    Spikes_test_centered = Spikes_test - Spikes_mean;

    % Train linear regression model on centered data
    mdl = fitlm(Spikes_train_centered', positions_train_centered);

    % Decode positions using the linear model on the centered test set
    decoded_positions_centered = predict(mdl, Spikes_test_centered');

    % Restore the decoded positions to the original scale by adding the mean of positions
    decoded_positions = decoded_positions_centered + positions_mean;

    % Calculate circular distance error
    num_positions = max(integer_pos);  % Assuming positions are from 1 to num_positions

    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), positions_test(i), num_positions);
        decoding_error = decoding_error + error^2;  % Note: Removed `round()` for better precision
    end

    % Mean squared circular error
    mean_squared_circular_error = decoding_error / length(decoded_positions);
end
