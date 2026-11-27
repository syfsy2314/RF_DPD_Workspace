clc;
clear;
close all;

% 滤波器带宽，单位 MHz
% BW = [40 50 60 70];
BW = [100 120 140 160 180];   % Nominal double-sided bandwidth, MHz
% H2 data order: [50 60 70 80 91] MHz (last point plotted at nominal 180 MHz)

% PA全带输出恢复信号的 NMSE，单位 dB
% 各方法的恢复 MAT 文件来自 test3 对应的 recovery 目录，模型使用 test2 的配对 IQ。
% 以 test2/y0.mat 为参考，对全部 38,400 个样点用 get_NMSE 复算。
% 汇总见 test3/result/test3_summary_report.md 的全带恢复 NMSE 章节。
NMSE_MG   = [-36.74 -37.48 -37.76 -37.80 -37.90];
NMSE_MR   = [-36.57 -37.50 -37.77 -38.24 -38.44];
NMSE_OPLS = [-37.86 -37.96 -38.41 -38.61 -38.64];
NMSE_SBL  = [-34.69 -35.93 -36.20 -36.58 -36.80];

figure('Color', 'w');

plot(BW, NMSE_MG,   '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
plot(BW, NMSE_MR,   '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(BW, NMSE_SBL,  '-d', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(BW, NMSE_OPLS, '-^', 'LineWidth', 1.8, 'MarkerSize', 7);

grid on;
box on;

xlabel('Acquisition Bandwidth (MHz)', 'FontSize', 13);
ylabel('NMSE(dB)', 'FontSize', 13);

legend({'MG', 'MR', 'SBL', 'LD-PLS'}, ...
       'Location', 'best', 'FontSize', 11);

% 横坐标从左到右递减
set(gca, 'XDir', 'reverse');
set(gca, 'XTick', BW);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

%xticks([50 60 70 80 90]);
% xlim([50 90]);
ylim([-39 -34]);

set(gca, 'FontSize', 12, 'LineWidth', 1.1, 'Color', 'w');

% title('不同带限滤波器带宽下PA全带输出恢复信号NMSE', 'FontSize', 13);
