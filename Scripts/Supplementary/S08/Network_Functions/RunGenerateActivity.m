
% Initialize variables and connections
memory_neurons_CA3 = false(1, N);
special_neurons_CA3=[];
MCells=find(memory_neurons_CA3==1);
recurrent_connections_CA3 = zeros(N, N);
IsTrue=zeros(N,T);

trajectory = levyFlightTrajectory(P, T);

time_bins_CA3 = zeros(TLearning, N);
firing_rate_snapshots_CA3 = cell(T, 1);
avg_activity_at_special_location_CA3 = zeros(N, 1);
Spikes = zeros(N,T);

special_radius_squared = special_radius^2;

for t = 1:T
    
    firing_rates_CA3 = SaveRates_CA3;
   
    % Apply recurrent connections to the firing rates
    updated_firing_rates_CA3 = firing_rates_CA3; 
    for i = 1:N
        if any(recurrent_connections_CA3(i, :) ~= 0)
            for j = 1:N
                if recurrent_connections_CA3(i, j) ~= 0
                    updated_firing_rates_CA3{i} = updated_firing_rates_CA3{i} + ...
                        recurrent_connections_CA3(i, j) * firing_rates_CA3{j};
                end
            end
        end
    end



% Add Gaussian noise to each 10x10 matrix independently
firing_rates_CA3 = cellfun(@(x) x + noise_std * randn(size(x)), updated_firing_rates_CA3, 'UniformOutput', false);    
firing_rates_CA3 = cellfun(@(x) max(x, 0), firing_rates_CA3, 'UniformOutput', false);


    peak_activities_CA3 = zeros(N, 1);

    % Find the peak activity 
    for i = 1:N
        firing_matrix = firing_rates_CA3{i};
        peak_activities_CA3(i) = max(firing_matrix(:));
    end

    normalized_firing_rates_CA3 = cell(size(firing_rates_CA3));
    for i = 1:N
        firing_matrix = firing_rates_CA3{i};
        if peak_activities_CA3(i) > 0
            normalized_firing_rates_CA3{i} = mu_rho * (firing_matrix / peak_activities_CA3(i));
        else
            normalized_firing_rates_CA3{i} = 0 .* firing_matrix;
        end
    end
    firing_rates_CA3 = normalized_firing_rates_CA3;



    % Get the current position 
    pos = trajectory(t, :);
    pos = round(pos); 

    active_neurons_CA3 = [];
    if all(pos > 0) && all(pos <= P) % Ensure the position is within bounds
        for i = 1:N
            if firing_rates_CA3{i}(pos(1), pos(2)) > thresholdCA3
                active_neurons_CA3 = [active_neurons_CA3, i];
            end
        end
    end

    spikes=find(cellfun(@(fr) fr(pos(1), pos(2)) > thresholdCA3, firing_rates_CA3));
    Spikes(spikes,t)=ones(length(spikes),1);

% Special location check 
if norm(pos - special_location) <= special_radius
    special_neurons_CA3 = [];

    % Count already identified memory cells
    num_memory_cells = sum(memory_neurons_CA3);

    for i = 1:N
        if num_memory_cells >= MemoryLimit
            break; 
        end

        if firing_rates_CA3{i}(pos(1), pos(2)) > mu_rho * threshold_goalCA3 && ~memory_neurons_CA3(i)
            special_neurons_CA3 = [special_neurons_CA3, i];
            memory_neurons_CA3(i) = true;
            num_memory_cells = num_memory_cells + 1; % Update count
        end
    end
end

SpecialCA3=find(memory_neurons_CA3);
if ~isempty(SpecialCA3)
    
    
    % Form strong connections among special neurons 
    for i = 1:length(SpecialCA3)
        A=SpecialCA3(i);
        for j = 1:length(SpecialCA3)
                    B=SpecialCA3(j);

          
                       if A ~= B

               
                recurrent_connections_CA3(A, B) = WeightCAMKII;
                IsTrue(A, t) = 1;

                        end

                        if A ~= B 
               
                recurrent_connections_CA3(A, B) = WeightCAMKII;
                IsTrue(B, t) = 1;

                        end

        end
    end
end

    if t > TLearning
        time_bins_CA3 = [time_bins_CA3(2:end, :); zeros(1, N)];
    end
    time_bins_CA3(mod(t - 1, TLearning) + 1, active_neurons_CA3) = 1;

    % Update recurrent connections 
    if t >= TLearning
        for i = active_neurons_CA3
            if time_bins_CA3(1, i) == 1 %
                for j = 2:TLearning
                    later_crossings = find(time_bins_CA3(j, :) == 1);
                    time_scale = (TLearning - j + 1) / TLearning; % Scale by temporal distance
                    for k = later_crossings
                        if i ~= k && (~memory_neurons_CA3(i) || ~memory_neurons_CA3(k)) 
                            recurrent_connections_CA3(i, k) = recurrent_connections_CA3(i, k) + WeightUpdatePlusCA3 * time_scale;
                            recurrent_connections_CA3(k, i) = recurrent_connections_CA3(k, i) - WeightUpdateMinusCA3 * time_scale; 
                        end
                    end
                end
            end
        end
    end


    % Apply decay to recurrent connections in CA3
    for i = 1:N
        for j = 1:N
           if ~memory_neurons_CA3(i) || ~memory_neurons_CA3(j) 
                recurrent_connections_CA3(i, j) = recurrent_connections_CA3(i, j) * (1 - decay_rate);
            end
        end
    end



% Normalize and threshold recurrent connections
max_val_CA3a = max(recurrent_connections_CA3(~memory_neurons_CA3, ~memory_neurons_CA3), [], 'all');
if max_val_CA3a > 0
    for i = 1:N
        for j = 1:N
            if ~memory_neurons_CA3(i) || ~memory_neurons_CA3(j) % Apply normalization only if not both are memory neurons
                recurrent_connections_CA3(i, j) = WeightScaleCA3 * (recurrent_connections_CA3(i, j) / max_val_CA3a);
            end
        end
    end
end

    recurrent_connections_CA3(isnan(recurrent_connections_CA3)) = 0;
    recurrent_connections_CA3((recurrent_connections_CA3./nanmax(recurrent_connections_CA3))<CA3Sparsity)=0;

    % Remove self-recurrent connections in CA3
    recurrent_connections_CA3(1:N + 1:end) = 0; 

    snapshot_CA3 = cell(N, 1);
    for i = 1:N
        snapshot_CA3{i} = firing_rates_CA3{i};
        if norm(pos - special_location) <= special_radius
            avg_activity_at_special_location_CA3(i) = avg_activity_at_special_location_CA3(i) + firing_rates_CA3{i}(special_location(1), special_location(2));
        end
    end
    firing_rate_snapshots_CA3{t} = snapshot_CA3;
end


avg_activity_at_special_location_CA3 = avg_activity_at_special_location_CA3 / T;

[~, sorted_indices_CA3] = sort(avg_activity_at_special_location_CA3, 'descend');


chosen_neurons_CA3 = sorted_indices_CA3(3:6);
chosen_neurons_CA3(3:4) = Choose2NeuronsCA3;


% Function to generate a PxP trajectory using Levy flight
function trajectory = levyFlightTrajectory(P, numTimePoints)
    trajectory = zeros(numTimePoints, 2);
    trajectory(1, :) = randi(P, 1, 2); % Starting position
    for t = 2:numTimePoints
        step_size = round(levyRandomStep());
        angle = 2 * pi * rand; % Random angle
        step = [step_size * cos(angle), step_size * sin(angle)];
        trajectory(t, :) = trajectory(t - 1, :) + step;
        % Ensure the position is within bounds
        trajectory(t, :) = max(min(trajectory(t, :), P), 1);
    end
end

% Function to generate a random step size using Levy distribution
function step = levyRandomStep()
    beta = 1.5; % Levy exponent
    scale_factor = 10; % Scaling factor to increase the step size
    u = randn * sqrt(gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)));
    v = randn;
    step = scale_factor * (u / abs(v)^(1 / beta));
end
