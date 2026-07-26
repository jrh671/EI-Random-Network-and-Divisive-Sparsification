function [NS2,ordered] = Compute_RateMap(Data,positions)

% Data=spike_mat_excit;
Figures = 0;

Spike = Data;

% Determining the path of the animal running through the track
integer_pos = floor(positions) + 1;         % the integer values of all the positions 

n_neu = size(Spike, 1);
    unique_positions = unique(integer_pos);
    n_pos = length(unique_positions);

    % Create a spatial matrix with rows for each neuron and columns for each position
    spatial_matrix = zeros(n_neu, n_pos);

    for p = 1:n_pos
        % Get the indices of the time bins where the position was 'p'
        pos_indices = integer_pos == unique_positions(p);

        % For 3
        % each neuron, sum up the spikes at these positions
        for i = 1:n_neu
            spatial_matrix(i, p) = mean(Spike(i, pos_indices));
        end
    end

   [~,PrefPos] = max(spatial_matrix,[],2);
    [~,ordered] = sort(PrefPos);

Ordering = ordered;
NS = spatial_matrix./max(abs(spatial_matrix));

NSt= spatial_matrix;
NS2 = exp(NSt)./max(exp(NSt),[],2);

figure;imagesc(NS2(Ordering,:),[0,1]);title('Position Ave Membrane Potential'); xlabel('Position');ylabel('Neuron (Sorted)'); colormap("hot");

colormap('hot')
colorbar



%% Normalize By Max Neuron Within Position
% RateMap=NS;

%% Normalize By Max Position Within Neurons
RateMap=NS2;

