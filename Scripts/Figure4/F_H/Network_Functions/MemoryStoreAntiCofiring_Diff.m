
for J=1:2

AAB=J-1;

if AAB==1
load('/Users/josehurtado/Documents/MATLAB/Manuscript/Home_Dir/Random_EIPlaceNet/Supplementary/SX/D_E_F/PreRun_Data/StoreMemory_AAB.mat')
else
load('/Users/josehurtado/Documents/MATLAB/Manuscript/Home_Dir/Random_EIPlaceNet/Supplementary/SX/D_E_F/PreRun_Data/StoreMemory_ABA.mat')
end

%% Compute Differences: (Memory - Non-Memory) within Each Trial, Context, and Epoch
numTrials = 55;
numEpochs = 6;

% Define target epochs
epoch4 = 4;
epoch6 = 6;

% Initialize accumulators for difference values
diff_context1_epoch4 = [];
diff_context2_epoch4 = [];
diff_context1_epoch6 = [];
diff_context2_epoch6 = [];

% Loop through trials
for Trial = 1:numTrials
    % Compute differences (Memory - Non-Memory) and store
    if epoch4 <= numEpochs
        diff_context1_epoch4 = [diff_context1_epoch4; StoreMemoryCells{Trial}(1, epoch4) - StoreNonMemoryCells{Trial}(1, epoch4)];
        diff_context2_epoch4 = [diff_context2_epoch4; StoreMemoryCells{Trial}(2, epoch4) - StoreNonMemoryCells{Trial}(2, epoch4)];
    end
    if epoch6 <= numEpochs
        diff_context1_epoch6 = [diff_context1_epoch6; StoreMemoryCells{Trial}(1, epoch6) - StoreNonMemoryCells{Trial}(1, epoch6)];
        diff_context2_epoch6 = [diff_context2_epoch6; StoreMemoryCells{Trial}(2, epoch6) - StoreNonMemoryCells{Trial}(2, epoch6)];
    end
end

%% Scatterplot: Epoch 6 vs. Epoch 8 (Difference Values with Connected Lines)
figure; hold on; 

% Define jitter amount
jitter_amount = 0; % Small jitter to prevent overlap

% Jittered x-coordinates
x1_c1 = diff_context1_epoch4 + jitter_amount * randn(size(diff_context1_epoch4));
y1_c1 = diff_context1_epoch6 + jitter_amount * randn(size(diff_context1_epoch6));

x2_c2 = diff_context2_epoch4 + jitter_amount * randn(size(diff_context2_epoch4));
y2_c2 = diff_context2_epoch6 + jitter_amount * randn(size(diff_context2_epoch6));

if AAB==1
% Scatter plot for differences in Context 1
scatter(x1_c1,y1_c1, 50, 'r', 'filled', 'MarkerFaceAlpha', 0.6); % Context 1
scatter(y2_c2,x2_c2, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6); % Context 2

else
% Scatter plot for differences in Context 1
scatter(y1_c1, x1_c1, 50, 'r', 'filled', 'MarkerFaceAlpha', 0.6); % Context 1
scatter(x2_c2, y2_c2, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6); % Context 2


end
% % Connect matching points within the same trial/context
% for i = 1:length(diff_context1_epoch4)
%     plot([diff_context1_epoch4(i), diff_context1_epoch6(i)], ...
%          [diff_context1_epoch4(i), diff_context1_epoch6(i)], 'r-', 'LineWidth', 1);
% end
% 
% for i = 1:length(diff_context2_epoch4)
%     plot([diff_context2_epoch4(i), diff_context2_epoch6(i)], ...
%          [diff_context2_epoch4(i), diff_context2_epoch6(i)], 'b-', 'LineWidth', 1);
% end

% Formatting
xlabel('First Epoch (Memory - Non-Memory)');
ylabel('Second Epoch (Memory - Non-Memory)');
if AAB==1
    title('Memory - Non-Memory (tA-AB & tB-BA)');
legend({'Context 2', 'Context 1'}, 'Location', 'Best');
else
    title('Memory - Non-Memory (tA-BA & tB-AB)');

legend({'Context 1', 'Context 2'}, 'Location', 'Best');

end
xlim([-0.1,1])
ylim([-0.1,1])
hold off;

%% Boxplots for Differences Across Conditions (WITH CONNECTED LINES)
figure; hold on;

% Boxplots: Independent for each condition
positions = [1, 2, 3, 4]; % X-axis positions
group_labels = {'Tag1 Neutral1', 'Tag1 Neutral2', 'Tag2 Neutral2', 'Tag2 Neutral1'};

boxplot(diff_context1_epoch4, 'positions', positions(1), 'Widths', 0.4, 'Symbol', 'o');
boxplot(diff_context1_epoch6, 'positions', positions(2), 'Widths', 0.4, 'Symbol', 'o');
boxplot(diff_context2_epoch4, 'positions', positions(3), 'Widths', 0.4, 'Symbol', 'o');
boxplot(diff_context2_epoch6, 'positions', positions(4), 'Widths', 0.4, 'Symbol', 'o');

jitter_amount = 0; % Small jitter to prevent overlap

% Overlay scatter points on top of boxplots
scatter_x1_c1 = positions(1) + jitter_amount * randn(size(diff_context1_epoch4));
scatter_x2_c1 = positions(2) + jitter_amount * randn(size(diff_context1_epoch6));
scatter_x1_c2 = positions(3) + jitter_amount * randn(size(diff_context2_epoch4));
scatter_x2_c2 = positions(4) + jitter_amount * randn(size(diff_context2_epoch6));

scatter(scatter_x1_c1, diff_context1_epoch4, 40, 'r', 'filled', 'MarkerFaceAlpha', 1);
scatter(scatter_x2_c1, diff_context1_epoch6, 40, 'r', 'filled', 'MarkerFaceAlpha', 1);
scatter(scatter_x1_c2, diff_context2_epoch4, 40, 'b', 'filled', 'MarkerFaceAlpha', 1);
scatter(scatter_x2_c2, diff_context2_epoch6, 40, 'b', 'filled', 'MarkerFaceAlpha', 1);


% Connect corresponding trials within the same context (with transparency)
for i = 1:length(diff_context1_epoch4)
    plot([scatter_x1_c1(i), scatter_x2_c1(i)], [diff_context1_epoch4(i), diff_context1_epoch6(i)], ...
         'Color', [1, 0, 0, 0.2], 'LineWidth', 0.5); % Red with 20% opacity
end

for i = 1:length(diff_context2_epoch4)
    plot([scatter_x1_c2(i), scatter_x2_c2(i)], [diff_context2_epoch4(i), diff_context2_epoch6(i)], ...
         'Color', [0, 0, 1, 0.2], 'LineWidth', 0.5); % Blue with 20% opacity
end

ylabel('Memory - NonMemory AntiCofiring Power')
if AAB==1
    title('Memory - Non-Memory (tA-AB & tB-BA)');
else
    title('Memory - Non-Memory (tA-BA & tB-AB)');


end

end
