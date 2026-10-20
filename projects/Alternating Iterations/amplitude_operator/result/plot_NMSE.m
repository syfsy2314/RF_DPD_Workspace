clc;
clear;
close all;

Rp_dB = [0.1 1 2 3 4 5 6];   %dB


NMSE_Ref   = [-39.11 -39.12 -39.08 -39.05 -39.05 -39.01 -39.01];
NMSE_Initial   = [-39.04 -37.92 -34.17 -31.36 -33.73 -27.15 -28.71];
NMSE_IterationFinal = [-38.54 -38.70 -37.47 -38.00 -36.60 -37.09 -36.88];

figure;

plot(Rp_dB, NMSE_Ref,   '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
plot(Rp_dB, NMSE_Initial,   '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
plot(Rp_dB, NMSE_IterationFinal,  '-d', 'LineWidth', 1.8, 'MarkerSize', 7);

grid on;
box on;

xlabel('Passband Ripple(dB)', 'FontSize', 13);
ylabel('NMSE(dB)', 'FontSize', 13);

legend({'Ref', 'Initial', 'IterationFinal', }, ...
       'Location', 'best', 'FontSize', 11);

% 横坐标从左到右递减
set(gca, 'XTick', Rp_dB);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

%xticks([50 60 70 80 90]);
% xlim([50 90]);
ylim([-40 -20]);

set(gca, 'FontSize', 12, 'LineWidth', 1.1);

% title('不同带限滤波器带宽下PA全带输出恢复信号NMSE', 'FontSize', 13);