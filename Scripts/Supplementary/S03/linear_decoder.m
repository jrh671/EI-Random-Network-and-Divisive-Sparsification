function [mean_squared_circular_error, decoded_positions] = linear_decoder(Spikes, integer_pos)

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
    mu = mean(Spikes_train, 2);  
    sigma = std(Spikes_train, 0, 2);  
    sigma(sigma == 0) = 1;  

    Spikes_train = (Spikes_train - mu) ./ sigma;  
    Spikes_test = (Spikes_test - mu) ./ sigma;  

    % Standardize and center the position labels
    mu_pos = mean(positions_train);
    sigma_pos = std(positions_train);
    if sigma_pos == 0
        sigma_pos = 1;  
    end
    positions_train = (positions_train - mu_pos) / sigma_pos;
    positions_test = (positions_test - mu_pos) / sigma_pos;

    % Train linear regression model
    mdl = fitlm(Spikes_train', positions_train);

    decoded_positions = predict(mdl, Spikes_test');

    % Unscale the decoded positions back to original scale
    decoded_positions = decoded_positions * sigma_pos + mu_pos;

    num_positions = max(integer_pos); 

    circular_distance = @(x, y, num_positions) min(abs(x - y), num_positions - abs(x - y));
    
    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), positions_test(i), num_positions);
        decoding_error = decoding_error + error^2;
    end
    
    % Mean squared circular error
    mean_squared_circular_error = decoding_error / length(decoded_positions);
end
