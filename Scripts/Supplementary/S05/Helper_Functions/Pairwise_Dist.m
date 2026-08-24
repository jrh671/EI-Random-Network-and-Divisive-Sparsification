%% ===== Pairwise distances on Figure 2 windows (FavPos; circular embedding) =====
% Requires: FavPos1, FavPos2, P  

assert(exist('FavPos1','var')==1 && exist('FavPos2','var')==1, 'FavPos1/2 not found.');
assert(exist('P','var')==1, 'P (position bins) not found.');

% Keep windows with no NaNs across neurons
valid1 = all(~isnan(FavPos1), 1);
valid2 = all(~isnan(FavPos2), 1);
assert(any(valid1) && any(valid2), 'No valid windows in one or both contexts.');

FP1 = FavPos1(:, valid1);   % [N x W1]
FP2 = FavPos2(:, valid2);   % [N x W2]

% --- Circular unit-circle embedding per neuron ---
favpos_embed = @(FP) [cos(2*pi*(FP-1)/P).',  sin(2*pi*(FP-1)/P).'];  % -> [W x 2N]

E1 = favpos_embed(FP1);   % [W1 x 2N]
E2 = favpos_embed(FP2);   % [W2 x 2N]

% --- Pairwise distances ---
d_within_1 = pdist(E1, 'euclidean');   
d_within_2 = pdist(E2, 'euclidean');   
D_WITHIN   = [d_within_1(:); d_within_2(:)];
D_ACROSS   = pdist2(E1, E2, 'euclidean');   
D_ACROSS   = D_ACROSS(:);



% --- Plot: two boxplots (Within vs Across) with colors, no grid ---
figure('Name','Figure 2: Pairwise distances (FavPos windows)','Color','w','Position',[140 140 720 460]);
h = boxplot([D_WITHIN; D_ACROSS], ...
            [repmat({'Within'}, numel(D_WITHIN), 1); repmat({'Across'}, numel(D_ACROSS), 1)], ...
            'Colors','k','Symbol','+'); 
ylabel('Euclidean distance (2N circular-embedded space)');
title(sprintf('Figure 2 windows (P=%d bins): Within-context vs Across-context', P));
grid off;

Sheet4a = D_WITHIN;

Sheet4b = D_ACROSS;

boxes = findobj(gca,'Tag','Box');
fillColors = [1 1 1; 0.5 0.5 0.5];  
for j = 1:numel(boxes)
    patch(get(boxes(j),'XData'), get(boxes(j),'YData'), ...
          fillColors(numel(boxes)-j+1,:), ...
          'FaceAlpha',0.5, 'EdgeColor','k', 'LineWidth',2);
end

ylim([0,35])

