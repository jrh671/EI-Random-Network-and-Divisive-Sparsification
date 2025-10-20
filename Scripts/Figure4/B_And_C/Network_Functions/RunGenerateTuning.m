if GenerateTuning==1

% Generate the Neuron vs. Input Matrix (A) for CA3
rng(0); %Ensure Same Connections Across Simulations
A_CA3 = alpha * rand(N, K);
rng(Seed); %Set Seed For Tuning Changes (Remapping)

% Initialize B as a cell array
B = cell(K, 1);
sigma = cell(N, 1); % Initialize sigma as a cell array
rho = cell(N, 1);   % Initialize rho as a cell array
Idx=[];
% Generate the Input vs. Position Matrix (B), sigma, and rho
for k = 1:K
    B{k} = zeros(P, P);
    indices = randperm(P * P, n_pos); % Randomly pick n_pos positions for tuning
    B{k}(indices) = 1;
    Idx(k,:)=indices;
end

SpecialPos=Idx;

% Compute the raw tuning matrix (M) for CA3
M_CA3 = cell(N, 1);
for i = 1:N
    M_CA3{i} = zeros(P, P);
    for j = 1:(P * P)
        [x, y] = ind2sub([P, P], j);
        for k = 1:K
            M_CA3{i}(x, y) = M_CA3{i}(x, y) + A_CA3(i, k) * B{k}(x, y);
        end
    end
end


% Normalize and threshold M for CA3
max_vals_CA3 = zeros(P, P);
for j = 1:(P * P)
    [x, y] = ind2sub([P, P], j);
    max_vals_CA3(x, y) = max(cellfun(@(M) M(x, y), M_CA3));
end

M_prim_CA3 = cell(N, 1);

for i = 1:N
sigma{i} = mu_sigma*ones(P, P); % Initialize with default sigma
rho{i} = mu_rho*ones(P, P);     % Initialize with default rho

    if sigma_sigma > 0
        if strcmp(distType, 'lognormal')
            sigma{i}(SpecialPos(i,:)) = lognrnd(mu_sigma, sigma_sigma, [1, n_pos]);
        else
            sigma{i}(SpecialPos(i,:)) = normrnd(mu_sigma, sigma_sigma, [1, n_pos]);
        end
    end

    if sigma_rho > 0
        if strcmp(distType, 'lognormal')
            rho{i}(SpecialPos(i,:)) = lognrnd(mu_rho, sigma_rho, [1, n_pos]);
        else
            rho{i}(SpecialPos(i,:)) = normrnd(mu_rho, sigma_rho, [1, n_pos]);
        end
    end

    M_prim_CA3{i} = zeros(P, P);
    for j = 1:(P * P)
        [x, y] = ind2sub([P, P], j);
        normalized_val = M_CA3{i}(x, y) / max_vals_CA3(x, y);
        if normalized_val >= theta
            M_prim_CA3{i}(x, y) = normalized_val;
        else
            M_prim_CA3{i}(x, y) = 0;
        end
    end
end


% Add noise and normalize again for CA3
M_prime_CA3 = cell(N, 1);
for i = 1:N
    Noise = alpha_p * pinknoise(P * P);
    Noise = reshape(Noise, [P, P]);
    M_prime_CA3{i} = M_prim_CA3{i} + Noise;
    M_prime_CA3{i} = M_prime_CA3{i} / max(M_prime_CA3{i}(:));
end

% Define the circular trajectory
[X, Y] = meshgrid(linspace(0, 360, P), linspace(0, 360, P));

% Compute the firing rates along the trajectory for CA3
firing_rates_CA3 = cell(N, 1);
for i = 1:N

    
    firing_rates_CA3{i} = zeros(P, P);
    for x = 1:P
        for y = 1:P
            for xi = 1:P
                for yi = 1:P
                    if M_prime_CA3{i}(xi, yi) > 0
                        if Periodic == 0
                            dist = (X(x, y) - X(xi, yi))^2 + (Y(x, y) - Y(xi, yi))^2;
                        else
                            dist = min(abs(X(x, y) - X(xi, yi)), 360 - abs(X(x, y) - X(xi, yi)))^2 + ...
                                   min(abs(Y(x, y) - Y(xi, yi)), 360 - abs(Y(x, y) - Y(xi, yi)))^2;
                        end
                        firing_rates_CA3{i}(x, y) = firing_rates_CA3{i}(x, y) + ...
                            rho{i}(xi, yi) * M_prime_CA3{i}(xi, yi) * exp(-dist / (2 * sigma{i}(xi, yi)^2));
                    end
                end
            end
        end
    end
end

% Initialize an array to store the peak activities of all neurons in CA3
peak_activities_CA3 = zeros(N, 1);

% Find the peak activity for each neuron in CA3
for i = 1:N
    firing_matrix = firing_rates_CA3{i};
    peak_activities_CA3(i) = max(firing_matrix(:));
end

% Find the global maximum peak activity across all neurons in CA3
global_peak_activity_CA3 = max(peak_activities_CA3);

normalized_firing_rates_CA3 = cell(size(firing_rates_CA3));
for i = 1:N
    firing_matrix = firing_rates_CA3{i};
    normalized_firing_rates_CA3{i} = mu_rho * (firing_matrix / global_peak_activity_CA3);
end

firing_rates_CA3 = normalized_firing_rates_CA3;
SaveRates_CA3 = firing_rates_CA3;

end