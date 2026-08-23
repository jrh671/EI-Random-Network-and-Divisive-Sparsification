for thresh = 1:length(thresholds)
    thresh
    for samples=1:N_samples
        
sigma = .8; % Standard deviation of the mixture of Gaussians  
Homogeneous = .05;
Rates = 4;
DivSparse = NormM;
DivSparse(DivSparse < thresholds(thresh)) = 0;

[NumNeurons, NumPositions] = size(DivSparse);
Pos = Positions;
Time = length(Pos);

Activity = zeros(NumNeurons, Time);


% Define homogeneous Poisson rates for each neuron
homogeneous_rates = rand(NumNeurons, 1) * Homogeneous; % Example: Random rates between 0 and 2

% Generate activity based on NormThresh and IntegerPos

for t = 1:Time
    posIdx = Pos(t);

    mixture_gaussians = zeros(NumNeurons, 1);

    for pos = 1:NumPositions

        neuron_amplitudes = DivSparse(:, pos);
        nonzero_indices = neuron_amplitudes > 0;

        % Create Gaussian envelopes for nonzero values
        for neuron = find(nonzero_indices)'
            amplitude = neuron_amplitudes(neuron);
            gaussian_envelope = amplitude * exp(-((1:NumPositions) - pos).^2 / (2 * sigma^2));
            mixture_gaussians(neuron) = mixture_gaussians(neuron) + gaussian_envelope(posIdx);
        end
    end

    % Add inhomogeneous and homogeneous Poisson activity
    Activity(:, t) = Rates.*poissrnd(mixture_gaussians) + poissrnd(homogeneous_rates);

end


X_DivSparse{thresh,samples} = Activity;
    end
end