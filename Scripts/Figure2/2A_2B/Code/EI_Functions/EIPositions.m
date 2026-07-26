%Setting up the pfs of the input neurons (continuous)
if NovelPF ==1
    rng(J) %First Arena
pf_cell1 = GET_TRACK_PFs(n_pos, prob_pf_center, n_input);
    rng((J+1e2)) %Second Arena
pf_cell2 = GET_TRACK_PFs(n_pos, prob_pf_center, n_input);
end


% Determining the speed of the animal running through the track
mean_dt = 0.005; sigma_dt = 0;              % time scale of updating activity
mean_v = 1; 
sigma_v = 0;      

dt_vec = normrnd(mean_dt, sigma_dt, 1, ceil(n_steps/mean_dt));
dt_vec(dt_vec < 0 ) = 0;
total_time_vec = cumsum(dt_vec);
v_vec0 = normrnd(mean_v, sigma_v, 1, ceil(n_steps/mean_dt)); 
v_vec = [0, v_vec0];
v_vec(v_vec < 0) = 0; 

% Determining the path of the animal running through the track
positions = 0.5 * (v_vec(1:end-1) + v_vec(2:end)) .* dt_vec; % generating the positions at each time step from dt and the difference in velocities (trying to make it smooth)
positions = mod(cumsum(positions), n_pos);
integer_pos = floor(positions) + 1;         % the integer values of all the positions 

