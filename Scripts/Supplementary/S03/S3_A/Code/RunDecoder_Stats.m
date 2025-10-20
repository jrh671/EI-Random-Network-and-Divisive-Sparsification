for idx = 1:(1 * 210)
    % Calculate P and N based on idx
    sample = mod(idx-1, 10) + 1;
    Ohm = ceil(idx / 10);

X{Ohm,sample}=results{idx};
end

excluded_p = 0;
time_bin_length=1;
n_pos=30;
n_laps=5;

positions=0:1:n_pos*n_laps-1;
for P=1:n_pos*n_laps
integer_pos(P)=mod(positions(P),n_pos)+1;
end

Length=1:21;


for L=1:10

    L

for N=Length

Spikes=X{N,L};

% 1: Plurality | 2: Template | 3: L1| 4: L2| 5: Linear
if Decoder == 1
    [decoding_error(N,L), neu_votes,decoded_pos] = plurality_voting_decoder(Spikes, integer_pos, excluded_p);
elseif Decoder == 2
    [decoding_error(N,L), decoded_positions, neu_votes] = template_matching_decoder(Spikes, integer_pos, excluded_p, time_bin_length);
elseif Decoder == 3
    [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l1(Spikes, integer_pos, .1);
elseif Decoder ==4
    [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, .1);
elseif Decoder ==5
    [decoding_error(N, L)] = linear_decoder(Spikes, integer_pos);
end



end
end

data = decoding_error; 

% Calculating the mean and standard error for each condition
means = mean(data, 2); % Mean across the second dimension
std_errors = std(data, [], 2) / sqrt(size(data, 2)); % Standard Error of the Mean

% Conditions (for x-axis)
if Instance<2
conditions = linspace(0,0.1,21); % Adjust this as per your conditions
elseif Instance>1
conditions = 0.0500*ones(1,21);
end

Colors = ['r','b','k','g','c'];

hold on;

    hErrorbar(Decoder) = errorbar(conditions(Length), means(Length), std_errors(Length), 'o', 'Color',  Colors(Decoder), 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');
   SaveMean{Decoder,1}= means;
    xlabel('Inhibitory Strength');
ylabel('Decoding Error');
title('Predicted vs Actual Position MSE');
ylim([0,100])