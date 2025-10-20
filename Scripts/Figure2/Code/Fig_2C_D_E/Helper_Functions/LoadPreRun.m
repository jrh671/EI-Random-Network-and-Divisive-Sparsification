load('Directory_To/PreRun_Data/resultsFast.mat')
load('Directory_To/PreRun_Data/W_InputEFast.mat')
load('Directory_To/PreRun_Data/PF_CellFast.mat')


pf_cell=pf_cellC;
W_inputE=W_inputECell{1};

n_pos=30;
n_laps=5;
mean_dt=.0005;
sigma_dt=0;
mean_v=1;
sigma_v=0;

%% Generate time and velocity vectors
n_steps = n_laps * n_pos;
dt_vec = normrnd(mean_dt, sigma_dt, 1, ceil(n_steps/mean_dt));
dt_vec(dt_vec < 0 ) = 0;
total_time_vec = cumsum(dt_vec);
v_vec = zeros(1,length(total_time_vec)+1);
v_vec0 = normrnd(mean_v, sigma_v, 1, ceil(n_steps/mean_dt)); 
v_vec(2:end) = v_vec0;
v_vec(v_vec < 0) = 0;

% Determining the path of the animal running through the track
positions = 0.5 * (v_vec(1:end-1) + v_vec(2:end)) .* dt_vec; % generating the positions at each time step from dt and the difference in velocities (trying to make it smooth)
positions = mod(cumsum(positions), n_pos);


N_params = 21; % Assuming 21 thresholds
N_samples = 10; % 10 samples per threshold

X1 = cell(N_params, N_samples); % Preallocate cell array

for i = 1:N_params
    for j = 1:N_samples
        X1{i, j} = results{(i-1)*N_samples + j}; % Assign each sample correctly
    end
end

X=X1;
pf_cell1=pf_cell;
Run_EffectiveTuning;
Indices1=neuron_indices;
Positions=ceil(SlidingAverage_s_mean(positions,mean_dt));

load('Directory_To/PreRun_Data/resultsSlow.mat')
load('Directory_To/PreRun_Data/W_InputESlow.mat')
load('Directory_To/PreRun_Data/PF_CellSlow.mat')

pf_cell=pf_cellC;
W_inputE=W_inputECell{1};

n_pos=30;
n_laps=5;
mean_dt=.005;
sigma_dt=0;
mean_v=1;
sigma_v=0;

%% Generate time and velocity vectors
n_steps = n_laps * n_pos;
dt_vec = normrnd(mean_dt, sigma_dt, 1, ceil(n_steps/mean_dt));
dt_vec(dt_vec < 0 ) = 0;
total_time_vec = cumsum(dt_vec);
v_vec = zeros(1,length(total_time_vec)+1);
v_vec0 = normrnd(mean_v, sigma_v, 1, ceil(n_steps/mean_dt)); 
v_vec(2:end) = v_vec0;
v_vec(v_vec < 0) = 0;

% Determining the path of the animal running through the track
positions = 0.5 * (v_vec(1:end-1) + v_vec(2:end)) .* dt_vec; % generating the positions at each time step from dt and the difference in velocities (trying to make it smooth)
positions = mod(cumsum(positions), n_pos);


N_params = 21; % Assuming 21 thresholds
N_samples = 10; % 10 samples per threshold

X2 = cell(N_params, N_samples); % Preallocate cell array

for i = 1:N_params
    for j = 1:N_samples
        X2{i, j} = results{(i-1)*N_samples + j}; % Assign each sample correctly
    end
end

X=X2;
pf_cell1=pf_cell;
Run_EffectiveTuning;
Indices2=neuron_indices;
Positions=ceil(SlidingAverage_s_mean(positions,mean_dt));