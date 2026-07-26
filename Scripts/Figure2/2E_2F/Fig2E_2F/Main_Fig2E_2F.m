%% ============================================================
% COMBINED FIGURE 3: Fast + Slow effective Vm compression
% Run from the directory where GainModulation.m is available
% ============================================================

addpath("AnalysisFunctions/")

Combine=0; %1 for Figure 2F. 0 For Individual

if Combine==0
Fast=1; %Fast Network (Set 0 for Slow time constant)

Load=1; %Load PreRun Simulation

if Load==1
if Fast==1
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2E_2F/PreRun_Data/Fast/InData.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2E_2F/PreRun_Data/Fast/OutData.mat')
else
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2E_2F/PreRun_Data/Slow/InData.mat')
load('/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2E_2F/PreRun_Data/Slow/OutData.mat')
end
end

InputSpikes=[InputSpikes1',InputSpikes2']';


GainModulation;
PlotGain;

else
   
base_path = '/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Internal/Figure2/2E_2F/PreRun_Data';

nbins = 20;

FastData = collect_effective_vm_data(fullfile(base_path,'Fast'));
SlowData = collect_effective_vm_data(fullfile(base_path,'Slow'));

%% -------------------------
% Common bin edges across both datasets
% -------------------------

X_all = [FastData.X; SlowData.X];
edges = linspace(min(X_all), max(X_all), nbins+1);

FastData = bin_effective_vm_data(FastData, edges);
SlowData = bin_effective_vm_data(SlowData, edges);

%% -------------------------
% Plot
% -------------------------

figure('Name','Effective Vm gain: Fast + Slow');
hold on

FastX=FastData.xc;
FastY=FastData.ym;
FastError=FastData.yse;
SlowX=SlowData.xc;
SlowY=SlowData.ym;
SlowError=SlowData.yse;

% Fast
errorbar(FastData.xc, FastData.ym, FastData.yse, ...
    'o-', ...
    'Color', [0 0.4470 0.7410], ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2, ...
    'CapSize', 8)

plot(FastData.xx, FastData.yy, ...
    '-', ...
    'Color', [0 0.4470 0.7410], ...
    'LineWidth', 3)

% Slow
errorbar(SlowData.xc, SlowData.ym, SlowData.yse, ...
    'o-', ...
    'Color', [0.8500 0.3250 0.0980], ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2, ...
    'CapSize', 8)

plot(SlowData.xx, SlowData.yy, ...
    '-', ...
    'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 3)

yline(0, 'k--', 'LineWidth', 1.5)
xline(0, 'k:', 'LineWidth', 1.5)

xlabel('Low-N Vm')
ylabel('High-N Vm - Low-N Vm')

title(sprintf(['Effective Vm compression: Fast \\alpha = %.3f, \\beta = %.3f; ', ...
               'Slow \\alpha = %.3f, \\beta = %.3f'], ...
    FastData.alpha_effective, FastData.beta_effective, ...
    SlowData.alpha_effective, SlowData.beta_effective))

box off

%% -------------------------
% Print summaries
% -------------------------

print_effective_vm_summary(FastData, 'Fast');
print_effective_vm_summary(SlowData, 'Slow');


end
%% ============================================================
% Local helper functions
% ============================================================

function Data = collect_effective_vm_data(data_dir)

    load(fullfile(data_dir,'InData.mat'))
    load(fullfile(data_dir,'OutData.mat'))

    InputSpikes = [InputSpikes1', InputSpikes2']';

    close all 
    GainModulation;

    valid_vm = results_Vm.included & ...
               isfinite(results_Vm.alpha) & ...
               isfinite(results_Vm.beta) & ...
               isfinite(results_Vm.R2_affine) & ...
               isfinite(results_Vm.rayleigh_R);

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

    Data.X = X;
    Data.Y = Y;
    Data.mdl = mdl;
    Data.m = m;
    Data.c = c;
    Data.alpha_effective = 1 + m;
    Data.beta_effective = c;
    Data.R2 = mdl.Rsquared.Ordinary;
    Data.N_neurons = sum(valid_vm);
    Data.N_points = numel(X);

end

function Data = bin_effective_vm_data(Data, edges)

    nbins = numel(edges)-1;

    xc = nan(nbins,1);
    ym = nan(nbins,1);
    yse = nan(nbins,1);
    n_per_bin = nan(nbins,1);

    for b = 1:nbins

        if b < nbins
            idx = Data.X >= edges(b) & Data.X < edges(b+1);
        else
            idx = Data.X >= edges(b) & Data.X <= edges(b+1);
        end

        n_per_bin(b) = sum(idx);

        if n_per_bin(b) < 5
            continue
        end

        xc(b) = mean(Data.X(idx), 'omitnan');
        ym(b) = mean(Data.Y(idx), 'omitnan');
        yse(b) = std(Data.Y(idx), 0, 'omitnan') ./ sqrt(n_per_bin(b));

    end

    xx = linspace(min(edges), max(edges), 200)';
    yy = predict(Data.mdl, table(xx, 'VariableNames', {'VmLow'}));

    Data.edges = edges;
    Data.xc = xc;
    Data.ym = ym;
    Data.yse = yse;
    Data.n_per_bin = n_per_bin;
    Data.xx = xx;
    Data.yy = yy;

end

function print_effective_vm_summary(Data, label)

    fprintf('\n===== %s effective Vm compression =====\n', label);
    fprintf('N neurons: %d\n', Data.N_neurons);
    fprintf('N pooled neuron-position points: %d\n', Data.N_points);
    fprintf('DeltaVm slope m: %.4f\n', Data.m);
    fprintf('DeltaVm intercept c: %.4f\n', Data.c);
    fprintf('Effective alpha = 1 + m: %.4f\n', Data.alpha_effective);
    fprintf('Effective beta = c: %.4f\n', Data.beta_effective);
    fprintf('R2: %.4f\n', Data.R2);

end