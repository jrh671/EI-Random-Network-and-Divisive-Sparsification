function mean_squared_circular_error = linear_decoder(Spikes, integer_pos)
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

    % Standardize and center training data
    mu = mean(Spikes_train, 2);  % Mean across training samples
    sigma = std(Spikes_train, 0, 2);  % Standard deviation across training samples
    sigma(sigma == 0) = 1;  % Avoid division by zero

    Spikes_train = (Spikes_train - mu) ./ sigma;  % Standardize training set
    Spikes_test = (Spikes_test - mu) ./ sigma;  % Apply same transformation to test set

    % Standardize and center the position labels
    mu_pos = mean(positions_train);
    sigma_pos = std(positions_train);
    if sigma_pos == 0
        sigma_pos = 1;  % Avoid division by zero
    end
    positions_train = (positions_train - mu_pos) / sigma_pos;
    positions_test = (positions_test - mu_pos) / sigma_pos;

    % Train linear regression model
    mdl = fitlm(Spikes_train', positions_train);

    % Decode positions using the linear model on the test set
    decoded_positions = predict(mdl, Spikes_test');

    % Unscale the decoded positions back to original scale
    decoded_positions = decoded_positions * sigma_pos + mu_pos;

    % Calculate circular distance error
    num_positions = max(integer_pos); % Assuming positions are from 1 to num_positions

    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), positions_test(i), num_positions);
        decoding_error = round(decoding_error + error^2);
    end

    % Mean squared circular error
    mean_squared_circular_error = decoding_error / length(decoded_positions);
end
