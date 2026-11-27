clc; clear; close all;

%% Bandwidth settings: normal order first
BW = [100 120 140 160 180];   % MHz

%% NMSE data, unit: dB
% 各方法曲线对应 test3 各方法 capture 目录下的实测 PA 输出回采 y_dpd_*.bin。
% 全带 MP 基准来自 test2/MP/y0_dpd.bin，NMSE 以对齐后的原始输入 x 为参考。
% 复算记录见 test3/result/test3_summary_report.md 的实测 DPD 章节。
% 当前绘图值中，MG 的 80、91 MHz 及 SBL 的 91 MHz 与复算记录有差异。
% H2 data order: [50 60 70 80 91] MHz (last point plotted at nominal 180 MHz)
NMSE_OPLS = [-37.12 -37.59 -37.76 -37.77 -37.86];
NMSE_MG   = [-35.29 -36.17 -36.33 -36.59 -36.81];
NMSE_MR   = [-35.76 -36.45 -36.82 -37.06 -37.24];
NMSE_SBL  = [-33.22 -34.52 -34.92 -35.27 -36.12];
NMSE_ref  = [-38.55 -38.55 -38.55 -38.55 -38.55];

%% Plot normally first
figure('Color', 'w');
%set(gcf, 'Position', [300 250 680 420]);

plot(BW, NMSE_MG,   '-o', 'LineWidth', 1.6, 'MarkerSize', 7); hold on;
plot(BW, NMSE_MR,   '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
plot(BW, NMSE_SBL,  '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
plot(BW, NMSE_OPLS, '-^', 'LineWidth', 1.6, 'MarkerSize', 7);
plot(BW, NMSE_ref, '-x', 'LineWidth', 1.6, 'MarkerSize', 7);

%% Axis and grid settings
ax = gca;
set(ax, 'Color', 'w');
set(ax, 'XDir', 'reverse');     % Reverse x-axis after normal plotting
set(ax, 'XTick', BW);
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12);

grid on;

xlabel('Acquisition Bandwidth (MHz)', 'FontName', 'Times New Roman', 'FontSize', 13);
ylabel('NMSE (dB)', 'FontName', 'Times New Roman', 'FontSize', 13);

legend('MG-DPD', 'MR-DPD', 'SBL-DPD', 'LD-PLS-DPD','FS-DPD Ref', ...
       'Location', 'best', ...
       'FontName', 'Times New Roman', ...
       'FontSize', 11);

ylim([-39 -31]);

%% Export high-resolution figure
% exportgraphics(gcf, 'NMSE_Trend.png', 'Resolution', 600);
