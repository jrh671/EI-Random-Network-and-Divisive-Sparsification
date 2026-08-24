%% Place-field peaks near special_location

nCompare = 20;
SPECIAL_RADIUS = 2;
JITTER_AMOUNT = 1;

MCells_use = MCells(1:min(nCompare, numel(MCells)));

all_neurons = 1:N;
eligible_random = setdiff(all_neurons, MCells);

if numel(eligible_random) < nCompare
    error('Not enough non-MCells available.');
end

rng('shuffle');
random_cells = eligible_random(randperm(numel(eligible_random), nCompare));

M_peak_at_special = false(numel(MCells_use), 1);
R_peak_at_special = false(numel(random_cells), 1);

M_peak_locations = nan(numel(MCells_use), 2);
R_peak_locations = nan(numel(random_cells), 2);

% MCells
for k = 1:numel(MCells_use)

    cell_id = MCells_use(k);
    map_idx = find(chosen_neurons_CA3ave == cell_id, 1);

    if isempty(map_idx)
        warning('MCell %d is not in chosen_neurons_CA3ave.', cell_id);
        continue
    end

    rate_map = cumulative_activity_map_CA3(:,:,map_idx);

    if all(~isfinite(rate_map(:))) || max(rate_map(:), [], 'omitnan') <= 0
        continue
    end

    [~, idx] = max(rate_map(:), [], 'omitnan');
    [peak_x, peak_y] = ind2sub(size(rate_map), idx);

    M_peak_locations(k,:) = [peak_x peak_y];

    d = sqrt((peak_x - special_location(1))^2 + ...
             (peak_y - special_location(2))^2);

    M_peak_at_special(k) = d <= SPECIAL_RADIUS;
end

% Random cells
for k = 1:numel(random_cells)

    cell_id = random_cells(k);
    map_idx = find(chosen_neurons_CA3ave == cell_id, 1);

    if isempty(map_idx)
        warning('Random cell %d is not in chosen_neurons_CA3ave.', cell_id);
        continue
    end

    rate_map = cumulative_activity_map_CA3(:,:,map_idx);

    if all(~isfinite(rate_map(:))) || max(rate_map(:), [], 'omitnan') <= 0
        continue
    end

    [~, idx] = max(rate_map(:), [], 'omitnan');
    [peak_x, peak_y] = ind2sub(size(rate_map), idx);

    R_peak_locations(k,:) = [peak_x peak_y];

    d = sqrt((peak_x - special_location(1))^2 + ...
             (peak_y - special_location(2))^2);

    R_peak_at_special(k) = d <= SPECIAL_RADIUS;
end

valid_M = all(isfinite(M_peak_locations), 2);
valid_R = all(isfinite(R_peak_locations), 2);

prop_M = mean(M_peak_at_special(valid_M));
prop_R = mean(R_peak_at_special(valid_R));

fprintf('\nWithin radius %.2f of special_location [%g %g]:\n', ...
    SPECIAL_RADIUS, special_location(1), special_location(2));

fprintf('MCells: %.3f (%d/%d)\n', ...
    prop_M, sum(M_peak_at_special(valid_M)), sum(valid_M));

fprintf('Random: %.3f (%d/%d)\n', ...
    prop_R, sum(R_peak_at_special(valid_R)), sum(valid_R));


%% Proportions
%% Plot

figure('Position', [200 200 1000 430]);

t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Proportions
nexttile;

bar([prop_R prop_M]);

set(gca, ...
    'XTick', 1:2, ...
    'XTickLabel', {'Random cells', 'Tagged'}, ...
    'FontSize', 11);

ylabel('Proportion of Cells (N = 20)');
ylim([0 1]);

title('Place Fields Near Goal Location');

box off;


% Peak locations
nexttile;
hold on;

rng(1);

nR = sum(valid_R);
nM = sum(valid_M);

R_jitter_x = (rand(nR,1) - 0.5) * 2 * JITTER_AMOUNT;
R_jitter_y = (rand(nR,1) - 0.5) * 2 * JITTER_AMOUNT;

M_jitter_x = (rand(nM,1) - 0.5) * 2 * JITTER_AMOUNT;
M_jitter_y = (rand(nM,1) - 0.5) * 2 * JITTER_AMOUNT;

scatter( ...
    R_peak_locations(valid_R,2) + R_jitter_x, ...
    R_peak_locations(valid_R,1) + R_jitter_y, ...
    60, 'o', 'filled');

scatter( ...
    M_peak_locations(valid_M,2) + M_jitter_x, ...
    M_peak_locations(valid_M,1) + M_jitter_y, ...
    80, '^', 'filled');

scatter( ...
    special_location(2), ...
    special_location(1), ...
    150, 'x', ...
    'LineWidth', 2);

theta = linspace(0, 2*pi, 300);

circle_x = special_location(2) + SPECIAL_RADIUS*cos(theta);
circle_y = special_location(1) + SPECIAL_RADIUS*sin(theta);

plot(circle_x, circle_y, 'k--', 'LineWidth', 1.5);

xlim([0 P]);
ylim([0 P]);

axis equal;
set(gca, 'YDir', 'normal', 'FontSize', 11);

xlabel('X Position');
ylabel('Y Position');

title('Place-field peak locations');

legend({'Random cells', 'Tagged', 'Special location', 'Radius'}, ...
    'Location', 'best');

box on;

sgtitle('Place-field peaks at special location', ...
    'FontWeight', 'bold');