function [DE, neu_votes, decoded_positions] = plurality_voting_decoder(Spikes, integer_pos, excluded_p, time_bin_length)


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

    % Calculate preferred firing location for each neuron using the first half data
    neu_votes = zeros(1, num_neurons);
    for i = 1:num_neurons
        if total_spikes(i) > threshold
            avg_fr_per_position = zeros(1, num_positions);
            for pos = 1:num_positions
                avg_fr_per_position(pos) = mean(Spikes_first_half(i, integer_pos_first_half == pos));
            end
            [~, neu_votes(i)] = max(avg_fr_per_position);
        end
    end

%     % Convert time bin length to number of 0.5 ms bins
%     time_bin_in_half_ms = 1;
% 
    % New number of time bins for the second half
    new_num_bins = length(Spikes_second_half(1,:));

%     % Sum the spikes for every time_bin_in_half_ms bins in the second half
%     new_Spikes = reshape(Spikes_second_half, [num_neurons, time_bin_in_half_ms, new_num_bins]);
%     new_Spikes = sum(new_Spikes, 2);
% 
%     % Take every time_bin_in_half_ms-th position in the second half
    new_integer_pos = integer_pos_second_half;

    % Initialize the decoded position vector
    decoded_positions = zeros(1, new_num_bins);

    % Go through each time bin
    for i = 1:new_num_bins
        % Initialize a votes vector to count votes for each position
        votes = zeros(1, num_positions);

        % Go through each neuron
        for j = 1:num_neurons
            % If the neuron spiked in this bin, has a preferred position and fires above the p-th percentile, add a vote for its preferred position
            if Spikes_second_half(j, i) > 0 && neu_votes(j) > 0
                votes(neu_votes(j)) = votes(neu_votes(j)) + 1;
            end
        end

        % The decoded position is the one with the most votes
        [~, decoded_positions(i)] = max(votes);
        decoded_positions(i)=round(decoded_positions(i));
    end


% Function to calculate circular distance
circular_distance = @(x, y) min(abs(x - y), num_positions - abs(x - y));

% Assuming 'decoded_positions' and 'new_integer_pos' are vectors of the same length
decoding_error = 0;
for i = 1:length(decoded_positions(1:length(Spikes_second_half(1,:))))
    error = circular_distance(decoded_positions(i), new_integer_pos(i));
    decoding_error = decoding_error + error^2;
end

DE= decoding_error/ length(decoded_positions);

   
