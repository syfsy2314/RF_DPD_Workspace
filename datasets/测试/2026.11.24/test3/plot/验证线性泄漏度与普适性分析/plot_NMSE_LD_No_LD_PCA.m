%% 验证线性泄漏与普适性分析
clc;
clear;
close all;

% 滤波器带宽，单位 MHz
BW = [100 120 140 160 180];   % MHz

% 原有三组数据，顺序对应 BW = [100 120 140 160 180]
NMSE_NESE_PCA  = [-33.57 -33.78 -34.76 -35.89 -36.34];


% LD-PLS 为 test3 归档的全带恢复结果，以 test2/y0.mat 为参考计算 NMSE。
% H2 单边带宽依次为 50、60、70、80、91 MHz，末点沿用标称双边 180 MHz。
NMSE_LD_PLS = [-37.86 -37.96 -38.41 -38.61 -38.64];

% 根据图片新增两组数据
NMSE_NESE_OPCA = [-35.90 -36.35 -37.02 -37.59 -37.73];

figure('Color', 'w');

plot(BW, NMSE_NESE_PCA, '-v', 'LineWidth', 1.8, 'MarkerSize', 7);hold;
plot(BW, NMSE_NESE_OPCA, '-d', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(BW, NMSE_LD_PLS, '-^', 'LineWidth', 1.8, 'MarkerSize', 7);

% 新增两项对比

grid on;
box on;

xlabel('Acquisition Bandwidth(MHz)', 'FontSize', 13);
ylabel('NMSE(dB)', 'FontSize', 13);

legend({'PCA','LD-PCA','LD-PLS' }, ...
       'Location', 'best', 'FontSize', 11);

% 横坐标从左到右递减
set(gca, 'XDir', 'reverse');

% 改这个可以改刻度
set(gca, 'XTick', BW);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'w');
ylim([-39 -32]);

% title('不同带限滤波器带宽下不同方法PA全带输出恢复信号NMSE', 'FontSize', 13);