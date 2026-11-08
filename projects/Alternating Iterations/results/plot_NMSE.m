clc;
clear;
close all;

% BW=80 新 H1 构造下的纹波扫描结果（标称纹波 0.1–1.9 dB，步进 0.2 dB）
Rp_dB = 0.1:0.2:1.9;   % dB

NMSE_Ref           = [-38.47 -38.47 -38.47 -38.46 -38.46 -38.46 -38.46 -38.45 -38.45 -38.45];
NMSE_Initial       = [-38.32 -37.86 -37.27 -36.60 -35.90 -35.19 -34.48 -33.78 -33.11 -32.45];
NMSE_IterationFinal = [-38.35 -38.28 -38.15 -37.97 -37.74 -37.47 -37.18 -36.86 -36.52 -36.18];

figure;

plot(Rp_dB, NMSE_Ref,            '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
plot(Rp_dB, NMSE_Initial,        '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Rp_dB, NMSE_IterationFinal, '-d', 'LineWidth', 1.8, 'MarkerSize', 7);

grid on;
box on;

xlabel('Passband Ripple (dB)', 'FontSize', 13);
ylabel('NMSE (dB)', 'FontSize', 13);

legend({'Ref', 'Initial', 'IterationFinal'}, ...
       'Location', 'best', 'FontSize', 11);

set(gca, 'XTick', Rp_dB);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1.1);

% 保存曲线到 results/ 目录
saveas(gcf, fullfile(fileparts(mfilename('fullpath')), 'NMSE_vs_ripple_BW80.png'));
