function [DE, decoded_positions, neu_votes] = template_matching_decoder(Spikes, integer_pos, excluded_p, time_bin_length)

    % Number of neurons
    num_neurons = size(Spikes, 1);

    % Number of time bins
    num_bins = size(Spikes, 2);

    % Number of positions
    num_positions = max(integer_pos);

    % Separate the Spikes and integer_pos into two halves
    Spikes_first_half = Spikes(:, 1:(num_bins/2));
    Spikes_second_half = Spikes(:, (num_bins/2 + 1):end);
    
    integer_pos_first_half = integer_pos(1:(num_bins/2));
    integer_pos_second_half = integer_pos((num_bins/2 + 1):end);

    % Calculate the total spike count for each neuron in the first half
    total_spikes = sum(Spikes_first_half, 2);

    % Find the firing rate threshold corresponding to the p-th percentile in the first half
    threshold = prctile(total_spikes, excluded_p);

    % Calculate average firing rate vectors for each position using the first half
    avg_fr_per_position = zeros(num_neurons, num_positions);
    for pos = 1:num_positions
        pos_indices = (integer_pos_first_half == pos);
        if any(pos_indices)
            avg_fr_per_position(:, pos) = mean(Spikes_first_half(:, pos_indices), 2);
        end
    end

    % Initialize the decoded position vector
    decoded_positions = zeros(1, length(Spikes_second_half(1,:)));

    % Go through each time bin in the second half
    for i = 1:length(Spikes_second_half(1,:))
        % Calculate the correlation between the current spike pattern and each average firing rate vector
        correlations = zeros(1, num_positions);
        for pos = 1:num_positions
            correlations(pos) = corr(Spikes_second_half(:, i), avg_fr_per_position(:, pos), 'Type', 'Pearson');
        end

        % The decoded position is the one with the highest correlation
        [~, decoded_positions(i)] = max(correlations);
    end

    % Calculate circular distance error
    circular_distance = @(x, y) min(abs(x - y), num_positions - abs(x - y));

    % Assuming 'decoded_positions' and 'integer_pos_second_half' are vectors of the same length
    decoding_error = 0;
    for i = 1:length(decoded_positions)
        error = circular_distance(decoded_positions(i), integer_pos_second_half(i));
        decoding_error = decoding_error + error^2;
    end

    % Mean squared circular error
    DE = decoding_error / length(decoded_positions);

    % Return neuron votes (most active position) for visualization or analysis
    neu_votes = avg_fr_per_position;
end
