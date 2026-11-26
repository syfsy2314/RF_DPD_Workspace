% 从工作区根目录运行 startup，然后运行本脚本。
% 依据 test2 的实测回采，计算 EVM 与左右邻道功率的线性平均。
scriptDir = fileparts(mfilename('fullpath'));
resultDir = fileparts(scriptDir);
testDir = fileparts(resultDir);
repoRoot = fileparts(fileparts(fileparts(fileparts(testDir))));
figureDir = fullfile(resultDir, 'figures');

bandwidthMHz = [40, 50, 60];
sampleCount = 38400;
[~, ~, x] = import_iq_binary( ...
    fullfile(repoRoot, 'datasets', 'data_38400_v0.6_test1.dat'), sampleCount);
x = x / max(abs(x));
evmValues = zeros(numel(bandwidthMHz), 3);
averageAcpr = zeros(size(evmValues));
for idx = 1:numel(bandwidthMHz)
    bw = bandwidthMHz(idx);
    captures = {fullfile(testDir, 'DOMP', sprintf('y_DOMP_bw%d_dpd.bin', bw)), ...
        fullfile(testDir, 'LD-PLS', sprintf('y_LD_PLS_bw%d_p_dpd.bin', bw)), ...
        fullfile(testDir, 'LD-PLS', sprintf('y_LD_PLS_bw%d_dpd.bin', bw))};
    for method = 1:numel(captures)
        [~, ~, y] = import_iq_binary(captures{method}, sampleCount);
        assert(numel(y) == sampleCount, '回采点数不匹配：%s', captures{method});
        [yAligned, xAligned] = align_circular_signal(y, x);
        evmValues(idx, method) = get_EVM(yAligned, xAligned);
        meter = comm.ACPR('SampleRate', 1e9, 'MainChannelFrequency', 0, ...
            'MainMeasurementBandwidth', 100e6, 'AdjacentMeasurementBandwidth', 100e6, ...
            'AdjacentChannelOffset', [-100e6, 100e6], ...
            'MainChannelPowerOutputPort', true, 'AdjacentChannelPowerOutputPort', true);
        acpr = step(meter, yAligned);
        averageAcpr(idx, method) = 10 * log10(mean(10.^(acpr / 10)));
        fprintf('BW %d method %d EVM %.9f ACPR %.9f %.9f AACPR %.9f\n', ...
            bw, method, evmValues(idx, method), acpr(1), acpr(2), averageAcpr(idx, method));
    end
end

plot_metric(bandwidthMHz, evmValues, 'test2 实测 DPD EVM', ...
    'EVM (%)', [1.20, 1.50], 1.20:0.05:1.50, ...
    fullfile(figureDir, 'test2_dpd_evm_vs_bandwidth.png'));
plot_metric(bandwidthMHz, averageAcpr, 'test2 实测 DPD 平均 ACPR', ...
    '平均 ACPR (dBc)', [-50, -45.5], -50:0.5:-45.5, ...
    fullfile(figureDir, 'test2_dpd_aacpr_vs_bandwidth.png'));

function plot_metric(bandwidthMHz, values, chartTitle, axisLabel, yRange, yTicks, outputPath)
fig = figure('Name', chartTitle, 'Color', 'white', 'Position', [100, 100, 880, 510]);
ax = axes(fig, 'FontName', 'Microsoft YaHei', 'FontSize', 12, ...
    'LineWidth', 0.9, 'YGrid', 'on', 'GridAlpha', 0.15, 'Box', 'off');
hold(ax, 'on');
linePenalty = plot(ax, bandwidthMHz, values(:, 2), '--s', ...
    'Color', [0.77, 0.31, 0.15], 'LineWidth', 1.9, 'MarkerSize', 10, ...
    'DisplayName', 'LD-PLS（有全带惩罚）');
lineDOMP = plot(ax, bandwidthMHz, values(:, 1), '-o', ...
    'Color', [0.12, 0.30, 0.55], 'LineWidth', 1.9, 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.12, 0.30, 0.55], 'DisplayName', 'DOMP');
lineNoPenalty = plot(ax, bandwidthMHz, values(:, 3), '-.^', ...
    'Color', [0.06, 0.49, 0.43], 'LineWidth', 1.9, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'white', 'DisplayName', 'LD-PLS（无全带惩罚）');
xticks(ax, bandwidthMHz);
xlim(ax, [37, 63]);
ylim(ax, yRange);
yticks(ax, yTicks);
xlabel(ax, 'H2 观测半带宽 (MHz)');
ylabel(ax, axisLabel);
title(ax, chartTitle, 'FontWeight', 'normal');
legend(ax, [lineDOMP, linePenalty, lineNoPenalty], 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
exportgraphics(fig, outputPath, 'BackgroundColor', 'white', 'Resolution', 300);
end
