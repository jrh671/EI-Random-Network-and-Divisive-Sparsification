nConditions = length(conditions);
nSamples = 10;

X = cell(nConditions,nSamples);
decoding_error = nan(nConditions,nSamples);

for idx = 1:(nConditions * nSamples)

    sample = mod(idx-1, nSamples) + 1;
    Ohm = ceil(idx / nSamples);

    X{Ohm,sample} = results{idx};

end

excluded_p = 0;
time_bin_length = 1;
n_pos = 30;
n_laps = 5;

positions = 0:1:n_pos*n_laps-1;
integer_pos = zeros(1,n_pos*n_laps);

for P = 1:n_pos*n_laps
    integer_pos(P) = mod(positions(P),n_pos) + 1;
end

Length = 1:nConditions;
lambda=0.1;

for L = 1:nSamples

    L

    for N = Length

        Spikes = X{N,L};

        % 1: Assembly Tagging | 2: Template | 3: L1 | 4: L2 | 5: Linear
        if Decoder == 1
            [decoding_error(N,L), neu_votes, decoded_pos] = ...
                plurality_voting_decoder(Spikes, integer_pos, excluded_p);

        elseif Decoder == 2
            [decoding_error(N,L), decoded_positions, neu_votes] = ...
                template_matching_decoder(Spikes, integer_pos, excluded_p, time_bin_length);

        elseif Decoder == 3
            [decoding_error(N,L), sparsity, prediction, TestPos] = ...
                linear_decoder_with_l1(Spikes, integer_pos, lambda);

        elseif Decoder == 4
            [decoding_error(N,L), sparsity, prediction, TestPos] = ...
                linear_decoder_with_l2(Spikes, integer_pos, lambda);

        elseif Decoder == 5
            decoding_error(N,L) = linear_decoder(Spikes, integer_pos);
        end

    end
end

data = decoding_error;

means = mean(data, 2, 'omitnan');
std_errors = std(data, [], 2, 'omitnan') ./ sqrt(sum(~isnan(data), 2));

Colors = ['r','b','k','g','c'];

hold on;

hErrorbar(Decoder) = errorbar(conditions, means, std_errors, ...
    'o', ...
    'Color', Colors(Decoder), ...
    'HandleVisibility', 'off');

SaveMean{Decoder,1} = means;

xlabel('Inhibitory Strength');
ylabel('Decoding Error');
title('Predicted vs Actual Position MSE');
ylim([0,100]);