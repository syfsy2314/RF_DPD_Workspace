clc;
clear;
close all;

% 从主实验保存的数组绘图，避免图和数值报告使用不同轮次的结果。
resultDir = fileparts(mfilename('fullpath'));
data = load(fullfile(resultDir, 'ripple_scan.mat'), ...
    'Rp_list', 'NMSE_Ref', 'NMSE_Init', 'NMSE_Iter');
Rp_dB = data.Rp_list;
NMSE_Ref = data.NMSE_Ref;
NMSE_Initial = data.NMSE_Init;
NMSE_IterationFinal = data.NMSE_Iter;

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

% 只导出坐标轴，避免把 MATLAB 图窗工具栏保存进结果图。
exportgraphics(gca, fullfile(resultDir, 'NMSE_vs_ripple_BW80.png'), ...
    'Resolution', 200);
