clc;
clear;
close all;

%% Data
% Acquisition bandwidth arranged from large to small
bandwidth = [180, 160, 140, 120, 100];

% Number of latent variables
ALVS_LD_PLS   = [10, 11, 11, 8, 8];
GS_Opt_LD_PLS = [12, 13, 12, 8, 6];

data = [ALVS_LD_PLS(:), GS_Opt_LD_PLS(:)];

%% Figure with MATLAB default window size
figure('Color', 'w');

%% Grouped bar chart
b = bar(data, 'grouped', 'BarWidth', 0.75);

% Blue and orange colors
b(1).FaceColor = [0.20, 0.45, 0.75];
b(2).FaceColor = [0.90, 0.45, 0.15];

b(1).EdgeColor = 'k';
b(2).EdgeColor = 'k';

b(1).LineWidth = 0.8;
b(2).LineWidth = 0.8;

hold on;

%% Numerical labels
for k = 1:numel(b)
    xPosition = b(k).XEndPoints;
    yPosition = b(k).YEndPoints;
    labels = string(b(k).YData);

    text(xPosition, yPosition + 0.15, labels, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10);
end

%% Axis labels
xlabel('Acquisition Bandwidth (MHz)', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12);

ylabel('Number of Latent Variables', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12);

%% X-axis tick labels
xticks(1:numel(bandwidth));
xticklabels(string(bandwidth));

%% Legend
legend({'ALVS-LD-PLS', 'GS-Opt-LD-PLS'}, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 10, ...
    'Box', 'off');

%% Axis formatting
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'TickDir', 'in', ...
    'Box', 'on', ...
    'Layer', 'top');

ylim([0, 14.5]);
yticks(0:2:14);

% Keep only horizontal grid lines
grid on;
set(gca, ...
    'XGrid', 'off', ...
    'YGrid', 'on', ...
    'GridLineStyle', ':', ...
    'GridAlpha', 0.25);
