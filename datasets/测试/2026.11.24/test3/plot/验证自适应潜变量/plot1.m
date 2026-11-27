clc;
clear;
close all;

%% ===================== 数据输入 =====================
% 反馈带宽，单位：MHz
% 数据顺序按照表格中的 90 MHz -> 50 MHz 排列
bandwidth = [100 120 140 160 180];   % MHz

% NMSE，单位：dB
ALVS   = fliplr([-37.5922 -37.9338 -37.6273 -37.0266 -36.5068]);
GS_Opt = fliplr([-37.9427 -38.1478 -37.7892 -37.0266 -36.8467]);
FL     = fliplr([-36.9855 -38.0133 -37.7703 -34.1623 -35.0803]);
FS     = fliplr([-34.9736 -34.9590 -34.8740 -34.6756 -34.7520]);

%% ===================== 绘图 =====================
figure;
hold on;
box on;
grid on;

plot(bandwidth, ALVS,   '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(bandwidth, GS_Opt, '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(bandwidth, FL,     '-^', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(bandwidth, FS,     '-d', 'LineWidth', 1.8, 'MarkerSize', 7);

%% ===================== 坐标轴设置 =====================
xlabel('Acquisition Bandwidth(MHz)', 'FontSize', 13);
ylabel('NMSE(dB)', 'FontSize', 13);

% 注意：xticks 必须是递增的数值向量
% xticks([50 60 70 80 90]);

% 设置坐标范围
% xlim([50 90]);

% 反转 x 轴，使其从左到右显示为 90 MHz -> 50 MHz
set(gca, 'XDir', 'reverse');
set(gca, 'XTick', bandwidth);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

ylim([-39 -32]);

%% ===================== 图例设置 =====================
legend({'ALVS-LD-PLS', ...
        'GS-Opt-LD-PLS', ...
        'FL-LD-PLS', ...
        'FS-LD-PLS'}, ...
        'Location', 'best', ...
        'FontSize', 11);

%% ===================== 图形美化 =====================
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.1);

% title('不同反馈带宽下恢复信号的NMSE对比', 'FontSize', 13);

hold off;
