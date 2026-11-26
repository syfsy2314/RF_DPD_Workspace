% test2 全带恢复 NMSE。数据来自相邻 recovery_nmse 目录的四份 test2_bw*_fullband_nmse.md。
bandwidthMHz = [30, 40, 50, 60];
nmseDOMP = [-37.194381321, -38.049079053, -37.877124893, -38.557225909];
nmseLDPLSPenalty = [-35.921252694, -37.836630373, -37.863192723, -37.964729741];
nmseLDPLSNoPenalty = [-35.494202897, -37.030985227, -36.944541393, -37.762853206];

fig = figure('Name', 'test2 全带恢复 NMSE', 'Color', 'white', ...
    'Position', [100, 100, 880, 510]);
ax = axes(fig, 'FontName', 'Microsoft YaHei', 'FontSize', 12, ...
    'LineWidth', 0.9, 'YGrid', 'on', 'GridAlpha', 0.15, 'Box', 'off');
hold(ax, 'on');

% 空心方形和实心圆使 50 MHz 处接近重合的两条曲线仍可区分。
linePenalty = plot(ax, bandwidthMHz, nmseLDPLSPenalty, '--s', ...
    'Color', [0.77, 0.31, 0.15], 'LineWidth', 1.9, 'MarkerSize', 10, ...
    'DisplayName', 'LD-PLS（有全带惩罚）');
lineDOMP = plot(ax, bandwidthMHz, nmseDOMP, '-o', ...
    'Color', [0.12, 0.30, 0.55], 'LineWidth', 1.9, 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.12, 0.30, 0.55], 'DisplayName', 'DOMP');
lineNoPenalty = plot(ax, bandwidthMHz, nmseLDPLSNoPenalty, '-.^', ...
    'Color', [0.06, 0.49, 0.43], 'LineWidth', 1.9, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'white', 'DisplayName', 'LD-PLS（无全带惩罚）');

xticks(ax, bandwidthMHz);
xlim(ax, [27, 63]);
ylim(ax, [-39, -35]);
yticks(ax, -39:0.5:-35);
xlabel(ax, 'H2 观测半带宽 (MHz)');
ylabel(ax, '全带恢复 NMSE (dB)');
title(ax, 'test2 不同带宽下的全带恢复 NMSE', 'FontWeight', 'normal');
legend(ax, [lineDOMP, linePenalty, lineNoPenalty], 'Location', 'southwest', ...
    'Box', 'off', 'FontSize', 11);

resultDir = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(resultDir, 'figures');
exportgraphics(fig, fullfile(outputDir, 'test2_recovery_nmse_vs_bandwidth.png'), ...
    'BackgroundColor', 'white', 'Resolution', 300);
