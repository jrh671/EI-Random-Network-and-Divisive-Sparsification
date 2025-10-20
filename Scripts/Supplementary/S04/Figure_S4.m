addpath("Helper_Functions/")

%% ====================== CONFIG & PATHS ======================
RATES      = 0;   % 1 = show frame-by-frame sorted imagesc during FavPos calc; 0 = skip
Save       = 0;   % 1 = save FavPos/RateMaps .mat files
PCA_DIMS   = 2;   % set to 2 or 3
LOAD_EXTRAS= 0;   % set 1 to load PF/Pos/W structs too

% NEW: choose which FavPos windows to animate (1-based). [] => plot all (original behavior).
PLOT_WINDOWS = [2 12 23];

p1_pf   = './SavedFiles/PF_CellFastSTDP1.mat';
p1_pos  = './SavedFiles/PosSTDP1.mat';
p1_res  = './SavedFiles/resultsFastSTDP1.mat';
p1_w    = './SavedFiles/W_InputEFastSTDP1.mat';

% Auto-derive STDP2 paths
p2_pf   = strrep(p1_pf,  'STDP1', 'STDP2');
p2_pos  = strrep(p1_pos, 'STDP1', 'STDP2');
p2_res  = strrep(p1_res, 'STDP1', 'STDP2');
p2_w    = strrep(p1_w,   'STDP1', 'STDP2');

%% ====================== LOAD ======================
if LOAD_EXTRAS
    S1_pf  = load(p1_pf);   %#ok<NASGU>
    S1_pos = load(p1_pos);  %#ok<NASGU>
    S1_w   = load(p1_w);    %#ok<NASGU>
    S2_pf  = load(p2_pf);   %#ok<NASGU>
    S2_pos = load(p2_pos);  %#ok<NASGU>
    S2_w   = load(p2_w);    %#ok<NASGU>
end

S1_res = load(p1_res);
S2_res = load(p2_res);

%% ====================== PARAMS ======================
P           = 30;   % positions per lap
window_size = 3;    % laps per window
step_size   = 2;    % advance by 2 laps

%% ====================== FAVPOS (STDP1) ======================
[~, X1_T] = size(S1_res.results);
assert(mod(X1_T, P)==0, 'STDP1 results: T must be a multiple of P (%d).', P);

[FavPos1, RateMaps1] = local_compute_FavPos_with_optional_rates( ...
    S1_res.results, P, window_size, step_size, RATES, ...
    'Sliding 2-lap windows (STDP1): neuron sorting by favorite position', ...
    PLOT_WINDOWS);

%% ====================== FAVPOS (STDP2) ======================
[~, X2_T] = size(S2_res.results);
assert(mod(X2_T, P)==0, 'STDP2 results: T must be a multiple of P (%d).', P);

[FavPos2, RateMaps2] = local_compute_FavPos_with_optional_rates( ...
    S2_res.results, P, window_size, step_size, RATES, ...
    'Sliding 2-lap windows (STDP2): neuron sorting by favorite position', ...
    PLOT_WINDOWS);

%% ===== (NEW) RATE MAPS: FIRST vs LAST THIRD (sorted by first valid STDP1 FavPos) =====
valid_col1 = find(all(~isnan(FavPos1),1), 1, 'first');
if ~isempty(valid_col1)
    % sort by STDP1 FavPos (first valid window)
    sort_key = FavPos1(:, valid_col1);
    [~, sort_idx_STDP1] = sort(sort_key, 'ascend');

    % ---- STDP1 ----
    X1 = S1_res.results;                     % [N x T1]
    L1 = size(X1,2) / P;                     % laps (safe because of assert above)
    n_third_1   = max(1, floor(L1/3));
    laps1_first = 1:n_third_1;
    laps1_last  = (L1 - n_third_1 + 1):L1;

    X1_3 = reshape(X1, size(X1,1), P, L1);   % [N x P x L1]
    R1_first = mean(X1_3(:,:,laps1_first), 3, 'omitnan');
    R1_last  = mean(X1_3(:,:,laps1_last ), 3, 'omitnan');
    cm = max(R1_first, [], 1); cm(cm==0)=1; R1_first = R1_first ./ cm;
    cm = max(R1_last , [], 1); cm(cm==0)=1; R1_last  = R1_last  ./ cm;

    % ---- STDP2 ----
    X2 = S2_res.results;                     % [N x T2]
    L2 = size(X2,2) / P;
    n_third_2   = max(1, floor(L2/3));
    laps2_first = 1:n_third_2;
    laps2_last  = (L2 - n_third_2 + 1):L2;

    X2_3 = reshape(X2, size(X2,1), P, L2);   % [N x P x L2]
    R2_first = mean(X2_3(:,:,laps2_first), 3, 'omitnan');
    R2_last  = mean(X2_3(:,:,laps2_last ), 3, 'omitnan');
    cm = max(R2_first, [], 1); cm(cm==0)=1; R2_first = R2_first ./ cm;
    cm = max(R2_last , [], 1); cm(cm==0)=1; R2_last  = R2_last  ./ cm;

    % ---- Plot 2x2 ----
    figRL = figure('Name','Rate maps — First vs Last Third (sorted by STDP1)', ...
                   'Color','w','Units','pixels','Position',[80 80 1200 820]);
    tl = tiledlayout(figRL, 2, 2, 'Padding','compact','TileSpacing','compact');

    nexttile(tl,1);
    imagesc(R1_first(sort_idx_STDP1,:), [0 1]);
    title(sprintf('STDP1 — First third (laps 1–%d of %d)', n_third_1, L1));
    xlabel('Position bin (1..P)'); ylabel('Neuron (sorted by STDP1 FavPos)'); colorbar; axis tight; colormap(parula);

    Sheet1a=R1_first(sort_idx_STDP1,:);

    nexttile(tl,2);
    imagesc(R1_last(sort_idx_STDP1,:), [0 1]);
    title(sprintf('STDP1 — Last third (laps %d–%d)', L1-n_third_1+1, L1));
    xlabel('Position bin (1..P)'); ylabel('Neuron (sorted by STDP1 FavPos)'); colorbar; axis tight; colormap(parula);
   
    Sheet1b=R1_last(sort_idx_STDP1,:);

    nexttile(tl,3);
    imagesc(R2_first(sort_idx_STDP1,:), [0 1]);
    title(sprintf('STDP2 — First third (laps 1–%d of %d)', n_third_2, L2));
    xlabel('Position bin (1..P)'); ylabel('Neuron (sorted by STDP1 FavPos)'); colorbar; axis tight; colormap(parula);

    Sheet1c=R2_first(sort_idx_STDP1,:);

    nexttile(tl,4);
    imagesc(R2_last(sort_idx_STDP1,:), [0 1]);
    title(sprintf('STDP2 — Last third (laps %d–%d)', L2-n_third_2+1, L2));
    xlabel('Position bin (1..P)'); ylabel('Neuron (sorted by STDP1 FavPos)'); colorbar; axis tight; colormap(parula);

    Sheet1d=R2_last(sort_idx_STDP1,:);

else
    warning('No valid FavPos1 window found; skipping first/last third rate-map plots.');
end
%% ===== END NEW BLOCK =====

if Save==1
    save(strrep(p1_res, 'resultsFastSTDP1.mat','FavPosFastSTDP1.mat'), 'FavPos1');
    save(strrep(p1_res, 'resultsFastSTDP1.mat','RateMapsFastSTDP1.mat'), 'RateMaps1');
    save(strrep(p2_res, 'resultsFastSTDP2.mat','FavPosFastSTDP2.mat'), 'FavPos2');
    save(strrep(p2_res, 'resultsFastSTDP2.mat','RateMapsFastSTDP2.mat'), 'RateMaps2');
end

%% ====================== JOINT PCA OF FavPos (2D/3D, no lines) ======================
valid1 = all(~isnan(FavPos1), 1);
valid2 = all(~isnan(FavPos2), 1);

F1 = FavPos1(:, valid1).';   % [nWin1 x N]
F2 = FavPos2(:, valid2).';   % [nWin2 x N]
if ~isempty(F1) && ~isempty(F2)
    assert(size(F1,2)==size(F2,2), 'FavPos feature dims mismatch between sessions.');
end

if size(F1,1) >= 2 && size(F2,1) >= 2
    Fjoint = [F1; F2];
    [~, score] = pca(Fjoint, 'Centered', true);
    Y = score(:, 1:max(3, PCA_DIMS));

    n1  = size(F1,1);
    n2  = size(F2,1);
    idx1 = 1:n1;
    idx2 = (n1+1):(n1+n2);

    lightBlue = [0.70, 0.85, 1.00];  darkBlue = [0.00, 0.10, 0.50];
    lightOrg  = [1.00, 0.85, 0.60];  darkOrg  = [0.70, 0.30, 0.00];
    t1 = linspace(0,1,n1)';  cmap1 = (1-t1).*lightBlue + t1.*darkBlue;
    t2 = linspace(0,1,n2)';  cmap2 = (1-t2).*lightOrg  + t2.*darkOrg;

    figure('Name', sprintf('Joint PCA (%dD) of FavPos: STDP1 vs STDP2', PCA_DIMS), 'Color','w');
    axes; hold on;

    if PCA_DIMS == 3
        h1 = scatter3(Y(idx1,1), Y(idx1,2), Y(idx1,3), 120, cmap1, 'filled', 'MarkerEdgeColor',[0.2 0.2 0.4]);
        h2 = scatter3(Y(idx2,1), Y(idx2,2), Y(idx2,3), 120, cmap2, 'filled', 'MarkerEdgeColor',[0.4 0.2 0.0]);
        xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
        grid off; axis vis3d; axis tight; view(35,25);
    else
        h1 = scatter(Y(idx1,1), Y(idx1,2), 120, cmap1, 'filled', 'MarkerEdgeColor',[0.2 0.2 0.4]);
        
        Sheet2a = [Y(idx1,1), Y(idx1,2)];

        h2 = scatter(Y(idx2,1), Y(idx2,2), 120, cmap2, 'filled', 'MarkerEdgeColor',[0.4 0.2 0.0]);
        
        Sheet2b = [Y(idx2,1), Y(idx2,2)];
        
        xlabel('PC 1'); ylabel('PC 2');
        grid off; axis tight;
    end
xlim([-100,100])
ylim([-120,100])

    title(sprintf('Joint PCA of FavPos (%dD): STDP1 (blue grad) vs STDP2 (orange grad)', PCA_DIMS));
    legend([h1, h2], {'STDP1','STDP2'}, 'Location','bestoutside');
else
    warning('Joint PCA of FavPos skipped: need at least 2 valid windows in both sessions.');
end

%% ===== MOMENT-TO-MOMENT PCA FROM RAW X (2D/3D; time-bin samples, no lines) =====
X1 = S1_res.results;  % [N x T1]
X2 = S2_res.results;  % [N x T2]
assert(size(X1,1)==size(X2,1), 'Neuron count mismatch between sessions (rows of X).');

X1_t = X1.';   % [T1 x N]
X2_t = X2.';   % [T2 x N]

X_joint_t = [X1_t; X2_t];     % [T1+T2 x N]
[~, score_t] = pca(X_joint_t, 'Centered', true);
Yt = score_t(:, 1:max(3, PCA_DIMS));
T1 = size(X1_t,1);  T2 = size(X2_t,1);
Yt1 = Yt(1:T1, :);
Yt2 = Yt(T1+1:T1+T2, :);

lightBlue = [0.70, 0.85, 1.00];  darkBlue = [0.00, 0.10, 0.50];
lightOrg  = [1.00, 0.85, 0.60];  darkOrg  = [0.70, 0.30, 0.00];
t1 = linspace(0,1,T1)';  cmap1_t = (1 - t1).*lightBlue + t1.*darkBlue;
t2 = linspace(0,1,T2)';  cmap2_t = (1 - t2).*lightOrg  + t2.*darkOrg;

figure('Name', sprintf('Moment-to-Moment PCA (%dD) of Raw X: STDP1 vs STDP2', PCA_DIMS), 'Color','w');
axes; hold on;

if PCA_DIMS == 3
    h1t = scatter3(Yt1(:,1), Yt1(:,2), Yt1(:,3), 100, cmap1_t, 'filled', 'MarkerEdgeColor',[0.2 0.2 0.4], 'MarkerFaceAlpha',0.9);
    h2t = scatter3(Yt2(:,1), Yt2(:,2), Yt2(:,3), 100, cmap2_t, 'filled', 'MarkerEdgeColor',[0.4 0.2 0.0], 'MarkerFaceAlpha',0.9);
    xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
    grid off; axis vis3d; axis tight; view(35,25);
else

    h1t = scatter(Yt1(:,1), Yt1(:,2), 100, cmap1_t, 'filled', 'MarkerEdgeColor',[0.2 0.2 0.4], 'MarkerFaceAlpha',0.9);
    
    Sheet3a = [Yt1(:,1), Yt1(:,2)];

    h2t = scatter(Yt2(:,1), Yt2(:,2), 100, cmap2_t, 'filled', 'MarkerEdgeColor',[0.4 0.2 0.0], 'MarkerFaceAlpha',0.9);
   
    Sheet3b = [Yt2(:,1), Yt2(:,2)];
    
    xlabel('PC 1'); ylabel('PC 2');
    grid off; axis tight;
end

title(sprintf('Moment-to-Moment PCA of Raw Activity (%dD): STDP1 (blue grad) vs STDP2 (orange grad)', PCA_DIMS));
legend([h1t, h2t], {'STDP1 time bins','STDP2 time bins'}, 'Location','bestoutside');


Pairwise_Dist;

%% ====================== LOCAL FUNCTION ======================
function [FavPos, RateMaps] = local_compute_FavPos_with_optional_rates(X, P, window_size, step_size, RATES, figName, plot_windows)
    % X: [N x T]; P positions/lap; RATES: 1 -> animate selected windows
    % plot_windows: [] => plot all windows (original behavior)

    [N, T] = size(X); %#ok<NASGU>
    assert(mod(T, P)==0, 'T must be a multiple of P.');
    L = T / P;

    starts = 1:step_size:(L - window_size + 1);
    nWin   = numel(starts);
    FavPos   = nan(N, nWin);
    RateMaps = nan(N, P, nWin);

    X3 = reshape(X, [N, P, L]);  % [neurons x positions x laps]

    if nargin < 7 || isempty(plot_windows)
        plot_idxs = 1:nWin;           % original behavior
    else
        plot_idxs = intersect(unique(plot_windows(:).'), 1:nWin);
    end

    make_fig = (RATES==1) && ~isempty(plot_idxs);
    if make_fig
        figure('Name', figName, 'Color','w');
    end

    for w = 1:nWin
        lap_start   = starts(w);
        chosen_laps = lap_start:(lap_start + window_size - 1);

        rate_map_sel = mean(X3(:,:,chosen_laps), 3, 'omitnan');  % [N x P]

        colmax = max(rate_map_sel, [], 1);
        colmax(colmax==0) = 1;
        rate_map_sel = rate_map_sel ./ colmax;

        RateMaps(:,:,w) = rate_map_sel;
        [~, fav_pos]    = max(rate_map_sel, [], 2);
        FavPos(:, w)    = fav_pos;

        if make_fig && ismember(w, plot_idxs)
            figure;
            [~, sort_by_fav_pos] = sort(fav_pos, 'ascend');
            clf;
            imagesc(X(sort_by_fav_pos, :), [0 1]);
            axis tight;
            colormap(parula);
            colorbar;
            xlabel('Position (1..30)');
            ylabel('Neuron (sorted by favorite position)');
            title(sprintf('%s — Laps %d–%d (window %d/%d)', figName, ...
                  chosen_laps(1), chosen_laps(end), w, nWin));
            set(gca,'XTick',1:P);
            drawnow;
            pause(0.6);
        end
    end
end
