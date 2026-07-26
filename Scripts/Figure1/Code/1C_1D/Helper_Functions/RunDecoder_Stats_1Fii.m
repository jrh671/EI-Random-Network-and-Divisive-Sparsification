% load('XtotT_array.mat')

if Instance==1
figure;
end
% hold on
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

if Instance<5
Length=1:21;
elseif Instance>5
Length=3:3:21;
% Length=1:21;
end

for L=1:10
    L

for N=Length
% % F=results{L+(N-1)*10};
% % clear FR
% % FR=SlidingAverage_s(F,.1);
% % FR=FR(:,1:150);
% % Spikes=FR;
Spikes=X{N,L};


% [decoding_error(N,L), decoded_positions, neu_votes] = bayesian_decoder(Spikes, integer_pos, excluded_p, time_bin_length);


% 1: Plurality | 2: Template | 3: L1| 4: L2| 5: Linear
if Decoder == 1
    [decoding_error(N,L), neu_votes,decoded_pos] = plurality_voting_decoder(Spikes, integer_pos, excluded_p, time_bin_length);
elseif Decoder == 2
    [decoding_error(N,L), decoded_positions, neu_votes] = template_matching_decoder(Spikes, integer_pos, excluded_p, time_bin_length);
elseif Decoder == 3
    [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l1(Spikes, integer_pos, 1);
elseif Decoder ==4
    [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, 1);
elseif Decoder ==5
    [decoding_error(N, L)] = linear_decoder2(Spikes, integer_pos);
end
% 
% if Decoder == 1
%     [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, 0);
% elseif Decoder == 2
%     [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, .1);
% elseif Decoder == 3
%     [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, 1);
% elseif Decoder ==4
%     [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, 5);
% elseif Decoder ==5
%     [decoding_error(N, L), sparsity,prediction,TestPos] = linear_decoder_with_l2(Spikes, integer_pos, 10);
% end

    % Instance=2;


% plot(integer_pos(16:90))
% hold on
% plot(decoded_pos,'LineWidth',3);
% hold off
% title(N*0.005)
% pause(2)
end
end
% 
% plot(linspace(0,.1,21),decoding_error/max(decoding_error),'LineWidth',2)
% % ylim([0 1])

% Example NxL matrix
% data = decoding_error./max(decoding_error); % 10 conditions (N=10) and 5 sets (L=5)
data = decoding_error; % 10 conditions (N=10) and 5 sets (L=5)

% Calculating the mean and standard error for each condition
means = mean(data, 2); % Mean across the second dimension
std_errors = std(data, [], 2) / sqrt(size(data, 2)); % Standard Error of the Mean

% Conditions (for x-axis)
if Instance<2
conditions = linspace(0,0.1,21); % Adjust this as per your conditions
elseif Instance>1
% conditions = 0.0500*ones(1,21);
conditions = linspace(0,0.1,21); % Adjust this as per your conditions
end

Colors = ['r','b','r','g','c'];

hold on;

if DecodePlot==1

%     hErrorbar = errorbar(conditions(Length), means(Length)./max(means), std_errors(Length)./max(means), '*', 'Color', colors(Instance, :), 'HandleVisibility', 'off');
    hErrorbar = errorbar(conditions(Length), means(Length)/max(means(Length)), std_errors(Length)/max(means(Length)), 'o', 'Color',  Colors(Instance), 'HandleVisibility', 'off');%colors(Instance, :), 'HandleVisibility', 'off');

%     
%     % Get the color of the errorbar
% %     color = get(hErrorbar, 'Color');
%     
%     % Plot a line with the same color
%     plot(conditions, means/max(means), ':', 'Color', colors, 'HandleVisibility', 'off');
   SaveMean{Instance,1}= means;
    xlabel('Inhibitory Strength');
ylabel('Decoding Error');
title('Predicted vs Actual Position MSE');
elseif DecodePlot==2
        hErrorbar = errorbar(conditions, gradient(means), std_errors, '*', 'Color', Colors(Decoder), 'HandleVisibility', 'off');
    
    % Get the color of the errorbar
%     color = get(hErrorbar, 'Color');
    
  %  Plot a line with the same color
%     plot(conditions, gradient(means), '--', 'Color', colors(Instance, :), 'HandleVisibility', 'off');
xlabel('Inhibitory Strength');
ylabel('Decoding Error');
title('Predicted vs Actual Position MSE');

elseif DecodePlot==0
Tuning_Sparsity_Stats
   
end


% % Add a legend entry for this plot
% 

% % Plotting
% % figure;
% for i = 1:10
% plot(conditions,decoding_error(:,i)/max(decoding_error(:,i)),'k.','MarkerSize',20);hold on
% end
% %errorbar(conditions, means, std_errors, 'o');
% xlabel('Inhibitory Strength');
% ylabel('Decoding Error');
% title('Predicted vs Acutal Position MSE');
% 
% % Customize as needed (colors, limits, etc.)
ylim([0,1])