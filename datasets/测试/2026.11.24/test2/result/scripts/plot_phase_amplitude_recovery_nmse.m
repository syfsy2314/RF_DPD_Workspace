% 三组均以全部 38,400 个样点计算全带 NMSE。
nmseDB = [-19.466444396, -20.207420106, -37.863783831];
caseLabels = {'开关 1 直接覆盖 · H2 ±80 MHz', ...
    '开关 2 幅相补偿 · H2 ±80 MHz', '开关 2 幅相补偿 · H2 ±44 MHz'};
colors = [0.41, 0.45, 0.52; 0.77, 0.31, 0.15; 0.77, 0.31, 0.15];

fig = figure('Name', 'test2 幅相补偿全带恢复 NMSE', ...
    'Color', 'white', 'Position', [100, 100, 880, 360]);
ax = axes(fig, 'FontName', 'Microsoft YaHei', 'FontSize', 12, ...
    'LineWidth', 0.9, 'XGrid', 'on', 'GridAlpha', 0.15, 'Box', 'off');
hold(ax, 'on');
for idx = 1:numel(nmseDB)
    row = numel(nmseDB) - idx + 1;
    scatter(ax, nmseDB(idx), row, 64, colors(idx, :), 'filled');
    text(ax, nmseDB(idx) + 0.58, row, sprintf('%.2f dB', nmseDB(idx)), ...
        'FontName', 'Microsoft YaHei', 'FontSize', 11, ...
        'Color', colors(idx, :), 'VerticalAlignment', 'middle');
end
yticks(ax, 1:numel(caseLabels));
yticklabels(ax, flip(caseLabels));
ylim(ax, [0.5, 3.5]);
xticks(ax, -40:5:-15);
xlim(ax, [-40, -15]);
xlabel(ax, '全带恢复 NMSE (dB)');
title(ax, 'test2 幅相补偿全带恢复 NMSE', 'FontWeight', 'normal');

resultDir = fileparts(fileparts(mfilename('fullpath')));
exportgraphics(fig, fullfile(resultDir, 'figures', 'test2_phase_amplitude_recovery_nmse.png'), ...
    'BackgroundColor', 'white', 'Resolution', 300);
