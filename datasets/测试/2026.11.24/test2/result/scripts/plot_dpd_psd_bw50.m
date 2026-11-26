% 从工作区根目录运行 startup，再运行本脚本。
% Welch 参数与 RMS 归一化沿用 projects/LD-PLS/main.m。
scriptDir = fileparts(mfilename('fullpath'));
resultDir = fileparts(scriptDir);
testDir = fileparts(resultDir);
datasetsDir = fileparts(fileparts(fileparts(testDir)));

fsMHz = 1000;
nfft = 1024;
window = hann(nfft);
sampleCount = 38400;
[~, ~, x] = import_iq_binary( ...
    fullfile(datasetsDir, 'data_38400_v0.6_test1.dat'), sampleCount);
x = x / max(abs(x));
captures = {fullfile(testDir, 'DOMP', 'y_DOMP_bw50_dpd.bin'), ...
    fullfile(testDir, 'LD-PLS', 'y_LD_PLS_bw50_p_dpd.bin'), ...
    fullfile(testDir, 'LD-PLS', 'y_LD_PLS_bw50_dpd.bin')};
signals = complex(zeros(sampleCount, numel(captures)));
for method = 1:numel(captures)
    [~, ~, y] = import_iq_binary(captures{method}, sampleCount);
    assert(numel(y) == sampleCount, '回采点数不匹配：%s', captures{method});
    signals(:, method) = align_circular_signal(y, x);
end
gainRms = rms(signals, 1) / rms(x);
signals = signals ./ gainRms;
[psd, frequencyMHz] = pwelch(signals, window, nfft / 2, nfft, fsMHz, 'centered');
psdDb = 10 * log10(100 * psd);
assert(isequal(size(psdDb), [nfft, 3]) && all(isfinite(psdDb), 'all'));

fig = figure('Name', 'test2 ±50 MHz 实测 DPD PSD', 'Color', 'white', ...
    'Position', [100, 100, 880, 510]);
ax = axes(fig, 'FontName', 'Microsoft YaHei', 'FontSize', 12, ...
    'LineWidth', 0.9, 'YGrid', 'on', 'GridAlpha', 0.15, 'Box', 'off');
hold(ax, 'on');
lineDOMP = plot(ax, frequencyMHz, psdDb(:, 1), '-', ...
    'Color', [0.12, 0.30, 0.55], 'LineWidth', 1.6, 'DisplayName', 'DOMP');
linePenalty = plot(ax, frequencyMHz, psdDb(:, 2), '--', ...
    'Color', [0.77, 0.31, 0.15], 'LineWidth', 1.6, 'DisplayName', 'LD-PLS（有全带惩罚）');
lineNoPenalty = plot(ax, frequencyMHz, psdDb(:, 3), '-.', ...
    'Color', [0.06, 0.49, 0.43], 'LineWidth', 1.6, 'DisplayName', 'LD-PLS（无全带惩罚）');
xlim(ax, [-500, 500]);
xticks(ax, -500:100:500);
xlabel(ax, '频率偏移 (MHz)');
ylabel(ax, '归一化功率谱密度 (dB)');
title(ax, 'test2 · H2 ±50 MHz 实测 DPD PSD', 'FontWeight', 'normal');
legend(ax, [lineDOMP, linePenalty, lineNoPenalty], 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
exportgraphics(fig, fullfile(resultDir, 'figures', 'test2_bw50_dpd_psd.png'), ...
    'BackgroundColor', 'white', 'Resolution', 300);
fprintf('Exported %d PSD traces on %d frequency bins\n', size(psd, 2), numel(frequencyMHz));
