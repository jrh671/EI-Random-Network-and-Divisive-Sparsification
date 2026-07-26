if Plasticity == 1
    TLearning = 15;  % Learning time bins
    thresholdCA3 = 0.35; % Activity threshold
    threshold_goalCA3 = 10000; 
    decay_rate = 0.35;   % Decay rate to bring weights back to zero

else
    TLearning = 1;  % Learning time bins
    thresholdCA3 = 10000; % Activity threshold
    threshold_goalCA3 = 10000;
    decay_rate = 1;   % Decay rate to bring weights back to zero
end

WeightScaleCA3 = 1;
WeightCAMKII = 2;
WeightUpdatePlusCA3 = .1;   % Weight scale for updates
WeightUpdateMinusCA3 = 1;   % Weight scale for updates
CA3Sparsity=0;

CLimit = 10;
CLimits=[CLimit,CLimit,CLimit,CLimit];

% Parameters for distribution
mu_sigma = 75;
sigma_sigma = 10; % Standard deviation (or log-std) for sigma
mu_rho = 10;     % Mean (or log-mean) for rho
sigma_rho = 0;   % Standard deviation (or log-std) for rho
alpha = 0.01;    % Scaling factor for the uniform distribution
Ratio = 0.2;     % Proportion of positions tuned
n_pos = floor((P * P) * Ratio); % Number of tuning centers per input
theta = 0.99999;    % Threshold for normalization
alpha_p = 0;     % Pink noise level
distType = 'normal'; % 'normal' or 'lognormal' for sigma and rho distributions
Periodic = 0;

