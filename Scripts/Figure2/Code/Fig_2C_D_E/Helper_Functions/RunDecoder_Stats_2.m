% Optional normalization
if NormThresh == 1
    X = X_NormThresh;
end

excluded_p       = 0;
time_bin_length  = 1;
n_pos            = 30;
n_laps           = 5;
Length           = 1:21;

% Decode across samples and conditions
for L = 1:N_samples
    for N = Length
        Spikes = X{N, L}(:, 1:end);
        switch Decoder
            case 1  % Plurality
                [decoding_error(N, L), neu_votes, decoded_pos] = ...
                    plurality_voting_decoder(Spikes, Positions, excluded_p);
            case 2  % Template
                [decoding_error(N, L), decoded_positions, neu_votes] = ...
                    template_matching_decoder(Spikes, Positions, excluded_p, time_bin_length);
            case 3  % L1
                [decoding_error(N, L), sparsity, prediction, TestPos] = ...
                    linear_decoder_with_l1(Spikes, Positions, 1);
            case 4  % L2
                [decoding_error(N, L), sparsity, prediction, TestPos] = ...
                    linear_decoder_with_l2(Spikes, Positions, 1);
            case 5  % Linear
                decoding_error(N, L) = linear_decoder(Spikes, Positions);
        end
    end
end

% Summary stats
data      = decoding_error;
means     = mean(data, 2);
std_errors = std(data, [], 2) / sqrt(size(data, 2));

SaveMeans{Overlap} = means;
SaveSEM{Overlap}   = std_errors;

% X-axis conditions
if Instance < 2
    conditions = linspace(0, 0.3, 21);
else
    conditions = 0.01500 * ones(1, 21);
end

Colors = ['b','r','k','g','c'];

if ShowSparsity == 1
    if NormThresh == 1
        xVals = thresholds(Length);
    else
        xVals = conditions(Length);
    end

    figure;
    if NormThresh == 1
        xSpars = thresholds(Length);
    else
        xSpars = conditions(Length);
    end
    plot(xSpars, Sparsity, 'o');
    if Overlap < 3
        xlabel('Inhibition');
    else
        xlabel('Threshold');
    end
    ylabel('Proportion Of Inactive Units On Average');
    title('Sparsity');

else
    % figure;
    hold on;
    errorbar(Sparsity(Length), means(Length), std_errors(Length), 'o-', ...
        'Color', Colors(Overlap), 'HandleVisibility', 'off');
    xlabel('Sparsity');
    ylabel('Decoding Error');
    title('Decoding Error vs Sparsity');
end
