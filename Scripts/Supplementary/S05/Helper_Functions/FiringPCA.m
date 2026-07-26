%% ====================== MOMENT-TO-MOMENT 3D PCA FROM RAW X ======================
% Observations = time bins, Features = neurons
X1 = S1_res.results;              % [N x T]
X2 = S2_res.results;              % [N x T]
[N1, T1] = size(X1); %#ok<NASGU>
[N2, T2] = size(X2); %#ok<NASGU>

% Time x Neurons for PCA (rows = samples)
X1_t = X1.';   % [T1 x N1]
X2_t = X2.';   % [T2 x N2]

% To make joint projection comparable, fit PCA on the concatenated time bins
X_joint_t = [X1_t; X2_t];  % [T1+T2 x N]  (assumes same N across sessions)
[coeff_t, score_t, ~] = pca(X_joint_t, 'Centered', true); %#ok<ASGLU>

Yt = score_t(:, 1:3);      % first 3 PCs of moment-to-moment activity
Yt1 = Yt(1:T1, :);
Yt2 = Yt(T1+1:T1+T2, :);

% Gradients per session (time goes light -> dark), reusing your palette idea
lightBlue = [0.70, 0.85, 1.00];  darkBlue = [0.00, 0.10, 0.50];  % STDP1
lightOrg  = [1.00, 0.85, 0.60];  darkOrg  = [0.70, 0.30, 0.00];  % STDP2
t1 = linspace(0,1,T1)';  cmap1_t = (1 - t1).*lightBlue + t1.*darkBlue;
t2 = linspace(0,1,T2)';  cmap2_t = (1 - t2).*lightOrg  + t2.*darkOrg;

figure('Name','Moment-to-Moment PCA (3D) of Raw X: STDP1 vs STDP2','Color','w');
ax = axes; hold(ax,'on');

% Use smaller markers to keep dense time points readable
h1t = scatter3(Yt1(:,1), Yt1(:,2), Yt1(:,3), 12, cmap1_t, 'filled', ...
               'MarkerEdgeColor',[0.2 0.2 0.4], 'MarkerFaceAlpha',0.9);
h2t = scatter3(Yt2(:,1), Yt2(:,2), Yt2(:,3), 12, cmap2_t, 'filled', ...
               'MarkerEdgeColor',[0.4 0.2 0.0], 'MarkerFaceAlpha',0.9);

xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
title('Moment-to-Moment PCA of Raw Activity (3D): STDP1 (blue grad) vs STDP2 (orange grad)');
grid off; axis vis3d; axis tight; view(35,25);
legend([h1t, h2t], {'STDP1 time bins','STDP2 time bins'}, 'Location','bestoutside');

% % Optional: save
% out_png = strrep(p1_res, 'resultsFastSTDP1.mat', 'MomentToMoment_PCA3D_STDP1_vs_STDP2.png');
% out_fig = strrep(p1_res, 'resultsFastSTDP1.mat', 'MomentToMoment_PCA3D_STDP1_vs_STDP2.fig');
% saveas(gcf, out_png); saveas(gcf, out_fig);
