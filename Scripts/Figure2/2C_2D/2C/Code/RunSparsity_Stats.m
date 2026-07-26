% Compute population sparsity as:
% proportion of neurons inactive at each time point,
% averaged across all time points

nConditions = length(conditions);
nSamples = 10;

X = cell(nConditions,nSamples);
sparsity_values = nan(nConditions,nSamples);

for idx = 1:(nConditions * nSamples)

    sample = mod(idx-1, nSamples) + 1;
    Ohm = ceil(idx / nSamples);

    X{Ohm,sample} = results{idx};

end

Length = 1:nConditions;

for L = 1:nSamples

    fprintf('Sample %d\n', L);

    for N = Length

        Spikes = X{N,L};

        % Assumes Spikes is neurons x time
        % sparsity at each time = fraction of neurons inactive
        sparsity_t = mean(Spikes == 0, 1);

        % average sparsity across all time bins
        sparsity_values(N,L) = mean(sparsity_t, 'omitnan');

        % If Spikes is time x neurons instead, use:
        % sparsity_t = mean(Spikes == 0, 2);
        % sparsity_values(N,L) = mean(sparsity_t, 'omitnan');

    end
end

data = sparsity_values;

means = mean(data, 2, 'omitnan');
std_errors = std(data, [], 2, 'omitnan') ./ ...
    sqrt(sum(~isnan(data), 2));

figure;
hold on;

errorbar(conditions, means, std_errors, ...
    'o-', ...
    'LineWidth', 2, ...
    'MarkerSize', 8);

xlabel('Inhibitory Strength');
ylabel('Mean Sparsity');
title('Population Sparsity vs Inhibitory Strength');
ylim([0 1]);
grid on;

SaveMean{1,1} = means;