%% 验证线性泄漏与普适性分析
clc;
clear;
close all;

% 滤波器带宽，单位 MHz
BW = [100 120 140 160 180];   % MHz

% 原有三组数据，顺序对应 BW = [100 120 140 160 180]
NMSE_PLS = fliplr([-37.08 -36.01 -35.21 -32.56 -32.05]);
NMSE_No_LD_PLS = fliplr([-37.84 -37.47 -36.70 -34.93 -34.15]);
% LD-PLS 为 test3 归档的全带恢复结果，以 test2/y0.mat 为参考计算 NMSE。
% H2 单边带宽依次为 50、60、70、80、91 MHz，末点沿用标称双边 180 MHz。
NMSE_LD_PLS = [-37.86 -37.96 -38.41 -38.61 -38.64];



figure('Color', 'w');

plot(BW, NMSE_PLS, '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
%plot(BW, NMSE_NESE_PCA, '-v', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(BW, NMSE_No_LD_PLS, '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
%plot(BW, NMSE_NESE_OPCA, '-d', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(BW, NMSE_LD_PLS, '-^', 'LineWidth', 1.8, 'MarkerSize', 7);

% 新增两项对比



grid on;
box on;

xlabel('Acquisition Bandwidth(MHz)', 'FontSize', 13);
ylabel('NMSE(dB)', 'FontSize', 13);

legend({'PLS', 'OR-PLS', 'LD-PLS'}, ...
       'Location', 'best', 'FontSize', 11);

% 横坐标从左到右递减
set(gca, 'XDir', 'reverse');

% 改这个可以改刻度
set(gca, 'XTick', BW);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'w');

ylim([-40 -30]);
% exportgraphics(gcf, 'NMSE_Trend1.jpg', 'Resolution', 300, 'BackgroundColor', 'white');
% title('不同带限滤波器带宽下不同方法PA全带输出恢复信号NMSE', 'FontSize', 13);
