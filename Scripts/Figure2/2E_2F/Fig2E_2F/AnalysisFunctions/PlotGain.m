%% ============================================================
% CORE Vm COMPRESSION FIGURES — MINIMAL VERSION
%
% Run after GainModulation with PositionSpecific = 1
%
% Requires:
%   Vm_low, Vm_high
%   results_Vm
%   n_pos_bins
% ============================================================

%% -------------------------
% Settings
% -------------------------
 
top_rayleigh_n = 150;

center_bin = ceil(n_pos_bins / 2);
x = 1:n_pos_bins;

valid_vm = results_Vm.included & ...
           isfinite(results_Vm.alpha) & ...
           isfinite(results_Vm.beta) & ...
           isfinite(results_Vm.R2_affine) & ...
           isfinite(results_Vm.rayleigh_R);

valid_ids = find(valid_vm);

[~, sort_idx] = sort(results_Vm.rayleigh_R(valid_vm), 'descend');
top_ids = valid_ids(sort_idx);
top_ids = top_ids(1:min(top_rayleigh_n, numel(top_ids)));

%% ============================================================
% FIGURE 1: Average aligned Vm tuning curves overlayed
% Gain-preserving normalization to Low-N amplitude
% ============================================================

aligned_low_vm = nan(numel(top_ids), n_pos_bins);
aligned_high_vm = nan(numel(top_ids), n_pos_bins);

for k = 1:numel(top_ids)

    n = top_ids(k);

    low = Vm_low(n,:);
    high = Vm_high(n,:);

    good = isfinite(low) & isfinite(high);

    if sum(good) < 3
        continue
    end

    [~, peak_bin] = max(low);

    shift_amt = center_bin - peak_bin;

    low_shift = circshift(low, shift_amt);
    high_shift = circshift(high, shift_amt);

    low_base = min(low_shift, [], 'omitnan');
    low_amp = max(low_shift - low_base, [], 'omitnan');

    if low_amp <= eps
        continue
    end

    aligned_low_vm(k,:) = (low_shift - low_base) / low_amp;
    aligned_high_vm(k,:) = (high_shift - low_base) / low_amp;

end

mean_low_vm = mean(aligned_low_vm, 1, 'omitnan');
sem_low_vm = std(aligned_low_vm, 0, 1, 'omitnan') ./ ...
             sqrt(sum(isfinite(aligned_low_vm), 1));

mean_high_vm = mean(aligned_high_vm, 1, 'omitnan');
sem_high_vm = std(aligned_high_vm, 0, 1, 'omitnan') ./ ...
              sqrt(sum(isfinite(aligned_high_vm), 1));

figure('Name','Average aligned Vm tuning curves');

fill([x fliplr(x)], ...
     [mean_low_vm + sem_low_vm, fliplr(mean_low_vm - sem_low_vm)], ...
     [0.85 0.85 0.85], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.35)

hold on

fill([x fliplr(x)], ...
     [mean_high_vm + sem_high_vm, fliplr(mean_high_vm - sem_high_vm)], ...
     [0.65 0.65 0.65], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.35)

plot(x, mean_low_vm, 'LineWidth', 3)
plot(x, mean_high_vm, 'LineWidth', 3)

xline(center_bin, 'k--', 'LineWidth', 1.5)
yline(1, 'k:', 'LineWidth', 1.5)

xlabel('Position bin aligned to Low-N peak')
ylabel('Vm normalized to Low-N amplitude')
title('Average aligned Vm tuning curves')
legend({'Low-N SEM','High-N SEM','Low-N','High-N'}, ...
       'Location','best')
ylim([0,1])
xlim([1,5])

%% ============================================================
% FIGURE 2: Aligned High-N minus Low-N difference
% Uses raw centered Vm values, not normalized values
% ============================================================

aligned_diff = nan(numel(top_ids), n_pos_bins);

for k = 1:numel(top_ids)

    n = top_ids(k);

    low = Vm_low(n,:);
    high = Vm_high(n,:);

    good = isfinite(low) & isfinite(high);

    if sum(good) < 3
        continue
    end

    [~, peak_bin] = max(low);

    shift_amt = center_bin - peak_bin;

    low_s = circshift(low, shift_amt);
    high_s = circshift(high, shift_amt);

    aligned_diff(k,:) = high_s - low_s;

end

mean_diff = mean(aligned_diff, 1, 'omitnan');
sem_diff = std(aligned_diff, 0, 1, 'omitnan') ./ ...
           sqrt(sum(isfinite(aligned_diff), 1));

figure('Name','Average aligned Vm difference');

fill([x fliplr(x)], ...
     [mean_diff + sem_diff, fliplr(mean_diff - sem_diff)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.4)

hold on

plot(x, mean_diff, 'LineWidth', 3)

yline(0, 'k--', 'LineWidth', 2)
xline(center_bin, 'k:', 'LineWidth', 2)

xlabel('Position bin aligned to Low-N peak')
ylabel('High-N Vm - Low-N Vm')
title('Average Vm suppression profile')

%% ============================================================
% FIGURE 3: Effective gain from DeltaVm vs VmLow
% ============================================================

X = [];
Y = [];

for n = find(valid_vm)'

    low  = Vm_low(n,:);
    high = Vm_high(n,:);

    good = isfinite(low) & isfinite(high);

    X = [X; low(good)'];
    Y = [Y; (high(good) - low(good))'];

end

tbl = table(X, Y, 'VariableNames', {'VmLow','DeltaVm'});

mdl = fitlm(tbl, 'DeltaVm ~ VmLow');

m = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, 'VmLow'));
c = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, '(Intercept)'));

m_p = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames, 'VmLow'));
c_p = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Properties.RowNames, '(Intercept)'));

alpha_effective = 1 + m;
beta_effective = c;

nbins = 20;
edges = linspace(min(X), max(X), nbins+1);

xc = nan(nbins,1);
ym = nan(nbins,1);
yse = nan(nbins,1);

for b = 1:nbins

    idx = X >= edges(b) & X < edges(b+1);

    if sum(idx) < 5
        continue
    end

    xc(b) = mean(X(idx), 'omitnan');
    ym(b) = mean(Y(idx), 'omitnan');
    yse(b) = std(Y(idx), 0, 'omitnan') ./ sqrt(sum(idx));

end

xx = linspace(min(X), max(X), 200)';
yy = predict(mdl, table(xx, 'VariableNames', {'VmLow'}));

figure('Name','Effective Vm gain');

errorbar(xc, ym, yse, 'o-', 'LineWidth', 2)
hold on

plot(xx, yy, 'LineWidth', 3)

yline(0, 'k--', 'LineWidth', 1.5)
xline(0, 'k:', 'LineWidth', 1.5)

xlabel('Low-N Vm')
ylabel('High-N Vm - Low-N Vm')
title(sprintf('Effective Vm compression: \\alpha = %.3f, \\beta = %.3f', ...
    alpha_effective, beta_effective))

legend({'Binned mean ± SEM', ...
        'Linear fit', ...
        'No change'}, ...
        'Location','best')

%% ============================================================
% Print compact summary
% ============================================================

fprintf('\n===== Core Vm compression summary =====\n');
fprintf('N neurons: %d\n', sum(valid_vm));
fprintf('Mean alpha from tuning fits: %.4f\n', ...
    mean(results_Vm.alpha(valid_vm), 'omitnan'));
fprintf('Median alpha from tuning fits: %.4f\n', ...
    median(results_Vm.alpha(valid_vm), 'omitnan'));
fprintf('Mean beta from tuning fits: %.4f\n', ...
    mean(results_Vm.beta(valid_vm), 'omitnan'));
fprintf('Median beta from tuning fits: %.4f\n', ...
    median(results_Vm.beta(valid_vm), 'omitnan'));
fprintf('Mean shape corr: %.4f\n', ...
    mean(results_Vm.shape_corr(valid_vm), 'omitnan'));
fprintf('Fraction alpha < 1: %.4f\n', ...
    mean(results_Vm.alpha(valid_vm) < 1, 'omitnan'));

fprintf('\nEffective gain from DeltaVm vs VmLow:\n');
fprintf('DeltaVm slope m: %.4f, p = %.3g\n', m, m_p);
fprintf('DeltaVm intercept c: %.4f, p = %.3g\n', c, c_p);
fprintf('Effective alpha = 1 + m: %.4f\n', alpha_effective);
fprintf('Effective beta = c: %.4f\n', beta_effective);
fprintf('R2: %.4f\n', mdl.Rsquared.Ordinary);