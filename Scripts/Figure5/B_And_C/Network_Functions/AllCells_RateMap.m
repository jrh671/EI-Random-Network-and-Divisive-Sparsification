% Identify neurons with nonzero tuning maps
zeroCells = cellfun(@(x) all(x(:) == 0), M_prim_CA3);
nonZeroIndices = find(~zeroCells);

RateMap = zeros(P, P, length(nonZeroIndices));

%% Compute occupancy

occupancy_map = zeros(P, P);

for t = 1:T

    pos = round(trajectory(t, :));

    if all(pos > 0) && all(pos <= P)
        occupancy_map(pos(1), pos(2)) = ...
            occupancy_map(pos(1), pos(2)) + 1;
    end

end

% Avoid division by zero in unvisited locations
occupancy_map(occupancy_map == 0) = 1;


%% Generate rate maps

for W = 1:length(nonZeroIndices)

    neuron_idx = nonZeroIndices(W);

    for t = 1:T

        pos = round(trajectory(t, :));

        if all(pos > 0) && all(pos <= P)

            firing_rates_CA3 = FiringRates{t}{neuron_idx};

            RateMap(pos(1), pos(2), W) = ...
                RateMap(pos(1), pos(2), W) + ...
                firing_rates_CA3(pos(1), pos(2));

        end

    end

    % Occupancy normalize 
    RateMap(:, :, W) = ...
        RateMap(:, :, W) ./ occupancy_map;

    % Normalize each neuron's map to its maximum
    max_val = max(RateMap(:, :, W), [], 'all');

    if max_val > 0
        RateMap(:, :, W) = ...
            RateMap(:, :, W) / max_val;
    end

end