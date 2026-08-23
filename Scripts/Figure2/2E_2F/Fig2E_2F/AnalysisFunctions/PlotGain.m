function Data = PlotGain(Vm_low, Vm_high, results_Vm, make_figures)
%PLOTGAIN Compute the effective-Vm analysis
%
% Data = PlotGain(Vm_low, Vm_high, results_Vm, make_figures)
%

if nargin < 4
    make_figures = true;
end

TOP_RAYLEIGH_N = 150;
NBINS = 20;

n_pos_bins = size(Vm_low, 2);
center_bin = ceil(n_pos_bins / 2);
x = 1:n_pos_bins;

valid_vm = results_Vm.included & ...
           isfinite(results_Vm.alpha) & ...
           isfinite(results_Vm.beta) & ...
           isfinite(results_Vm.R2_affine) & ...
           isfinite(results_Vm.rayleigh_R);

%% Effective Vm regression data
X = [];
Y = [];

for n = find(valid_vm)'
    low = Vm_low(n,:);
    high = Vm_high(n,:);
    good = isfinite(low) & isfinite(high);

    X = [X; low(good)']; %#ok<AGROW>
    Y = [Y; (high(good) - low(good))']; %#ok<AGROW>
end

if isempty(X)
    error('No valid Vm observations are available for the effective-gain analysis.');
end

tbl = table(X, Y, 'VariableNames', {'VmLow','DeltaVm'});
mdl = fitlm(tbl, 'DeltaVm ~ VmLow');

m = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, 'VmLow'));
c = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, '(Intercept)'));

Data.X = X;
Data.Y = Y;
Data.mdl = mdl;
Data.alpha_effective = 1 + m;
Data.beta_effective = c;
Data.valid_vm = valid_vm;

if ~make_figures
    return
end

%% Figure 1: average aligned Vm tuning curves
valid_ids = find(valid_vm);
[~, sort_idx] = sort(results_Vm.rayleigh_R(valid_vm), 'descend');
top_ids = valid_ids(sort_idx);
top_ids = top_ids(1:min(TOP_RAYLEIGH_N, numel(top_ids)));

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
     [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.35)
hold on
fill([x fliplr(x)], ...
     [mean_high_vm + sem_high_vm, fliplr(mean_high_vm - sem_high_vm)], ...
     [0.65 0.65 0.65], 'EdgeColor','none', 'FaceAlpha',0.35)
plot(x, mean_low_vm, 'LineWidth', 3)
plot(x, mean_high_vm, 'LineWidth', 3)
xline(center_bin, 'k--', 'LineWidth', 1.5)
yline(1, 'k:', 'LineWidth', 1.5)
xlabel('Position bin aligned to Low-N peak')
ylabel('Vm normalized to Low-N amplitude')
title('Average aligned Vm tuning curves')
legend({'Low-N SEM','High-N SEM','Low-N','High-N'}, 'Location','best')
ylim([0,1])
xlim([1,n_pos_bins])

%% Figure 2: effective gain from DeltaVm vs VmLow
Data = bin_effective_vm_data(Data, NBINS, []);

figure('Name','Effective Vm gain');
errorbar(Data.xc, Data.ym, Data.yse, 'o-', 'LineWidth', 2)
hold on
plot(Data.xx, Data.yy, 'LineWidth', 3)
yline(0, 'k--', 'LineWidth', 1.5)
xline(0, 'k:', 'LineWidth', 1.5)
xlabel('Low-N Vm')
ylabel('High-N Vm - Low-N Vm')
title(sprintf('Effective Vm compression: \\alpha = %.3f, \\beta = %.3f', ...
    Data.alpha_effective, Data.beta_effective))
legend({'Binned mean ± SEM','Linear fit','No change'}, 'Location','best')
end

function Data = bin_effective_vm_data(Data, nbins, edges)

if isempty(edges)
    edges = linspace(min(Data.X), max(Data.X), nbins + 1);
else
    nbins = numel(edges) - 1;
end

xc = nan(nbins,1);
ym = nan(nbins,1);
yse = nan(nbins,1);

for b = 1:nbins
    idx = Data.X >= edges(b) & Data.X < edges(b+1);

    n = sum(idx);
    if n < 5
        continue
    end

    xc(b) = mean(Data.X(idx), 'omitnan');
    ym(b) = mean(Data.Y(idx), 'omitnan');
    yse(b) = std(Data.Y(idx), 0, 'omitnan') / sqrt(n);
end

xx = linspace(min(edges), max(edges), 200)';
yy = predict(Data.mdl, table(xx, 'VariableNames', {'VmLow'}));

Data.edges = edges;
Data.xc = xc;
Data.ym = ym;
Data.yse = yse;
Data.xx = xx;
Data.yy = yy;
end
