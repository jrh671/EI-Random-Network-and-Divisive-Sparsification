%% ============================================================
% FIGURE 3: EFFECTIVE Vm COMPRESSION
%
%   - GainModulation.m
%   - PlotGain.m contains the shared effective-gain calculation
%   - CompareLowHigh.m plots the low and high input drive states
%
% Set COMBINE to false for one condition or true for Fast + Slow overlay.
% ============================================================

COMBINE = true;
FAST = true;  % used only when COMBINE is false

base_path = '/Users/josehurtado/Documents/MATLAB/Final_Manuscript/Random_EIPlaceNet_Data/Figure2/2E_2F/PreRun_Data';

if ~COMBINE
    if FAST
        data_dir = fullfile(base_path, 'Fast');
    else
        data_dir = fullfile(base_path, 'Slow');
    end

    load(fullfile(data_dir, 'InData.mat'))
    load(fullfile(data_dir, 'OutData.mat'))
    InputSpikes = [InputSpikes1', InputSpikes2']';

    GainModulation;
    PlotGain(Vm_low, Vm_high, results_Vm, true);

else
    FastData = collect_effective_vm_data(fullfile(base_path, 'Fast'));
    SlowData = collect_effective_vm_data(fullfile(base_path, 'Slow'));

    % Use identical display bins for both conditions so Fast and Slow
    % are compared over exactly the same Low-N Vm ranges.
    NBINS = 20;
    X_all = [FastData.X; SlowData.X];
    edges = linspace(min(X_all), max(X_all), NBINS + 1);

    FastData = bin_for_display(FastData, edges);
    SlowData = bin_for_display(SlowData, edges);

    % Compute x-axis composition separately for Fast and Slow.
    FastComposition = compute_xaxis_composition(FastData, edges);
    SlowComposition = compute_xaxis_composition(SlowData, edges);

    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    %% ============================================================
    % Combined effective Vm compression
    % ============================================================

    figure('Name','Effective Vm gain: Fast + Slow');
    hold on

    errorbar(FastData.xc, FastData.ym, FastData.yse, ...
        'o-', 'Color',[0 0.4470 0.7410], ...
        'MarkerFaceColor','none', 'LineWidth',2, 'CapSize',8)
    plot(FastData.xx, FastData.yy, '-', ...
        'Color',[0 0.4470 0.7410], 'LineWidth',3)

    errorbar(SlowData.xc, SlowData.ym, SlowData.yse, ...
        'o-', 'Color',[0.8500 0.3250 0.0980], ...
        'MarkerFaceColor','none', 'LineWidth',2, 'CapSize',8)
    plot(SlowData.xx, SlowData.yy, '-', ...
        'Color',[0.8500 0.3250 0.0980], 'LineWidth',3)

    yline(0,'k--','LineWidth',1.5)
    xline(0,'k:','LineWidth',1.5)
    xlabel('Low-N Vm')
    ylabel('High-N Vm - Low-N Vm')
    title(sprintf(['Effective Vm compression: Fast \\alpha = %.3f, \\beta = %.3f; ', ...
                   'Slow \\alpha = %.3f, \\beta = %.3f'], ...
        FastData.alpha_effective, FastData.beta_effective, ...
        SlowData.alpha_effective, SlowData.beta_effective))
    box off

    %% ============================================================
    % Low-N Vm x-axis composition: Fast and Slow separately
    % ============================================================

    figure('Name','Low-N Vm x-axis composition: Fast and Slow');

    tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    nexttile
    bar(bin_centers, ...
        [FastComposition.proportion_neurons, ...
         FastComposition.proportion_neuron_positions], ...
        'grouped');
    ylabel('Proportion (%)')
    title('Fast')
    legend({'Neurons with \geq1 position in bin (%)', ...
            'All neuron-position observations (%)'}, ...
            'Location','best', ...
            'Interpreter','tex')
    xlim([edges(1), edges(end)])
    box off

    nexttile
    bar(bin_centers, ...
        [SlowComposition.proportion_neurons, ...
         SlowComposition.proportion_neuron_positions], ...
        'grouped');
    xlabel('Low-N Vm')
    ylabel('Proportion (%)')
    title('Slow')
    xlim([edges(1), edges(end)])
    box off
end

if COMBINE == false
    CompareLowHigh
end

%% ============================================================
% Local helpers
% ============================================================

function Data = collect_effective_vm_data(data_dir)
    load(fullfile(data_dir, 'InData.mat'))
    load(fullfile(data_dir, 'OutData.mat'))
    InputSpikes = [InputSpikes1', InputSpikes2']';

    GainModulation;
    Data = PlotGain(Vm_low, Vm_high, results_Vm, false);

    % Keep neuron and position-group identities for each pooled Vm
    % observation.
    valid_vm = Data.valid_vm;
    NeuronID = [];
    PositionID = [];

    for n = find(valid_vm)'
        good = isfinite(Vm_low(n,:)) & isfinite(Vm_high(n,:));
        positions = find(good);

        NeuronID = [NeuronID; repmat(n, numel(positions), 1)]; %#ok<AGROW>
        PositionID = [PositionID; positions(:)]; %#ok<AGROW>
    end

    if numel(NeuronID) ~= numel(Data.X)
        error('Observation identities do not match the effective-gain data.');
    end

    Data.NeuronID = NeuronID;
    Data.PositionID = PositionID;
end

function Composition = compute_xaxis_composition(Data, edges)
    nbins = numel(edges) - 1;

    total_neurons = numel(unique(Data.NeuronID));
    total_neuron_positions = numel(Data.X);

    proportion_neurons = nan(nbins,1);
    proportion_neuron_positions = nan(nbins,1);

    for b = 1:nbins
        if b < nbins
            idx = Data.X >= edges(b) & Data.X < edges(b+1);
        else
            idx = Data.X >= edges(b) & Data.X <= edges(b+1);
        end

        proportion_neurons(b) = ...
            100 * numel(unique(Data.NeuronID(idx))) / total_neurons;

        proportion_neuron_positions(b) = ...
            100 * sum(idx) / total_neuron_positions;
    end

    Composition.proportion_neurons = proportion_neurons;
    Composition.proportion_neuron_positions = proportion_neuron_positions;
end

function Data = bin_for_display(Data, edges)
    nbins = numel(edges) - 1;
    xc = nan(nbins,1);
    ym = nan(nbins,1);
    yse = nan(nbins,1);

    for b = 1:nbins
        if b < nbins
            idx = Data.X >= edges(b) & Data.X < edges(b+1);
        else
            idx = Data.X >= edges(b) & Data.X <= edges(b+1);
        end

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

    Data.xc = xc;
    Data.ym = ym;
    Data.yse = yse;
    Data.xx = xx;
    Data.yy = yy;
end
