clc; clear; close all;

%% Bandwidth settings: normal order first
BW = [100 120 140 160 180];   % MHz

%% Average ACPR data, unit: dBc
% 各方法曲线对应 test3 各方法 capture 目录下的实测 PA 输出回采 y_dpd_*.bin。
% test2/MP/y0_dpd.bin 是全带 MP 基准回采，复算的平均 ACPR 为 -49.59 dBc。
% 复算记录见 test3/result/test3_summary_report.md 的实测 DPD 章节。
% 当前绘图值中，LD-PLS 的 50、70 MHz，MG 的 80 MHz 及 MP 基准与复算记录有差异。
% H2 data order: [50 60 70 80 91] MHz (last point plotted at nominal 180 MHz)
% 10*log10(mean(10.^(ACPR/10), 2)), rounded after averaging
ACLR_LDPLS = [-47.77 -48.02 -49.07 -49.27 -49.57];
ACLR_MG    = [-44.40 -46.46 -46.97 -47.28 -47.44];
ACLR_MR    = [-45.01 -46.78 -48.40 -48.64 -49.17];
ACLR_SBL   = [-39.63 -41.79 -42.94 -44.09 -44.47];
ACLR_ref   = [-50.53 -50.53 -50.53 -50.53 -50.53];

%% Plot ACLR trend
figure('Color', 'w');
% set(gcf, 'Position', [300 250 680 420]);

plot(BW, ACLR_MG,    '-o', 'LineWidth', 1.6, 'MarkerSize', 7); hold on;
plot(BW, ACLR_MR,    '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
plot(BW, ACLR_SBL,   '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
plot(BW, ACLR_LDPLS, '-^', 'LineWidth', 1.6, 'MarkerSize', 7);
plot(BW, ACLR_ref, '-x', 'LineWidth', 1.6, 'MarkerSize', 7);

%% Axis and grid settings
ax = gca;
set(ax, 'Color', 'w');
set(ax, 'XDir', 'reverse');     % Reverse x-axis after normal plotting
set(ax, 'XTick', BW);
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12);

grid on;


xlabel('Acquisition Bandwidth (MHz)', 'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('ACPR (dBc)', 'FontName', 'Times New Roman', 'FontSize', 13);

legend('MG-DPD', 'MR-DPD', 'SBL-DPD', 'LD-PLS-DPD','FS-DPD Ref',...
       'Location', 'best', ...
       'FontName', 'Times New Roman', ...
       'FontSize', 11);

ylim([-52 -35]);

%% Export high-resolution figure
% exportgraphics(gcf, 'ACLR_Trend.png', 'Resolution', 600);
