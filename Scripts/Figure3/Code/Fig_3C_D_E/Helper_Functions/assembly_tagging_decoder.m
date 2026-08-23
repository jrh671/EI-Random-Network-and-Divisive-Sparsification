function [DE, neu_votes, decoded_positions] = assembly_tagging_decoder(Spikes, Positions, excluded_p)

    
    num_neurons = size(Spikes, 1);

    num_bins = size(Spikes, 2);

    num_positions = max(Positions);

    % Separate the Spikes and integer_pos into two halves
    Spikes_first_half = Spikes(:, 1:(num_bins/2));
    Spikes_second_half = Spikes(:, (num_bins/2 + 1):end);
    
    integer_pos_first_half = Positions(1:(num_bins/2));
    integer_pos_second_half = Positions((num_bins/2 + 1):end);

    total_spikes = sum(Spikes_first_half, 2);
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


    new_num_bins = length(Spikes_second_half(1,:));

    new_integer_pos = integer_pos_second_half;

    decoded_positions = zeros(1, new_num_bins);

    for i = 1:new_num_bins
        votes = zeros(1, num_positions);

        for j = 1:num_neurons

            if Spikes_second_half(j, i) > 0 && neu_votes(j) > 0
                votes(neu_votes(j)) = votes(neu_votes(j)) + 1;
            end
        end

        % The decoded position is the one with the most votes
        [~, decoded_positions(i)] = max(votes);
        decoded_positions(i)=round(decoded_positions(i));
    end


circular_distance = @(x, y) min(abs(x - y), num_positions - abs(x - y));

decoding_error = 0;
for i = 1:length(decoded_positions(1:length(Spikes_second_half(1,:))))
    error = circular_distance(decoded_positions(i), new_integer_pos(i));
    decoding_error = decoding_error + error^2;
end

DE= decoding_error/ length(decoded_positions);

   
