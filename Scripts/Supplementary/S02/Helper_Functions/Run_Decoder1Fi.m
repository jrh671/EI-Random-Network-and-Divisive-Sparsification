
excluded_p = 0;

Spikes=FiringRate;

Decoder=1;

Positions = ceil(SlidingAverage_s_mean(integer_pos,mean_dt));

% 1: Plurality | 2: Template | 3: L1| 4: L2| 5: Linear
if Decoder == 1
    [decoding_error, neu_votes,decoded_pos] = plurality_voting_decoder(Spikes, Positions, excluded_p);
elseif Decoder == 2
    [decoding_error, decoded_positions, neu_votes] = template_matching_decoder(Spikes, Positions, excluded_p);
elseif Decoder == 3
    [decoding_error, sparsity,prediction,TestPos] = linear_decoder_with_l1(Spikes, Positions, 1);
elseif Decoder ==4
    [decoding_error, sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, Positions, 1);
elseif Decoder ==5
    [decoding_error] = linear_decoder2(Spikes, integer_pos);
end

data(Iteration,Option) = decoding_error; 

