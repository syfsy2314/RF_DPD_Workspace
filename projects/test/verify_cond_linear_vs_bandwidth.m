%% 猜想验证：带宽降低 -> 条件数升高（共线性病态）+ 线性分量占比上升
% 条件数用列归一化口径 kappa(H_n)，纯反映共线性，排除列尺度塌缩的干扰
clear; close all;
fs = 1000;
%% 导入数据（只加载一次）
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
%% 扫描设置
BW_list = 100:-5:30;             % 采集带宽（MHz），从大到小
nBW = numel(BW_list);
kappa_n = zeros(1, nBW);         % 列归一化条件数（共线性口径）
p_lin   = zeros(1, nBW);         % 线性分量功率占比
K = 7;
M = 5;
Dstop1 = 0.0001;
Dstop2 = 0.000001;
dens  = 20;
Dpass1 = (10^(0.5/20) - 1) / (10^(0.5/20) + 1);      % H1 通带波纹 0.5 dB
Dpass2 = 0.0057501127785;                            % H2 通带波纹
%% 带宽扫描
for iBW = 1:nBW
    BW = BW_list(iBW);
    % 观测链路滤波器（通带随 BW 缩放）
    [N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(fs/2), [1 0], [Dpass1, Dstop1]);
    H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
    [N2, Fo, Ao, W] = firpmord([BW, BW+3]/(fs/2), [1 0], [Dpass2, Dstop2]);
    H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
    h2 = H2_FIR.Numerator;
    % 观测链路
    y_f = H1_FIR(y0);
    [y_f, x_al, ~] = align_circular_signal(y_f, x);
    y_bf = H2_FIR(y_f);
    % 列归一化条件数（纯共线性）
    H_x = geth_even_BLMP(x_al, K, M, h2);
    Hn = H_x ./ vecnorm(H_x);
    s = svd(Hn, 'econ');
    kappa_n(iBW) = s(1) / s(end);
    % 线性分量功率占比
    x_lin = H_x(:, 1);
    w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);
    y_lin = x_lin * w_lin;
    p_lin(iBW) = (y_lin' * y_lin) / (y_bf' * y_bf);
    fprintf('BW = %3d MHz: cond(H_n) = %.3g, 线性占比 = %.2f%%\n', ...
        BW, kappa_n(iBW), 100 * p_lin(iBW));
end
%% 双纵轴：条件数与线性占比（横轴带宽从大到小）
figure;
yyaxis left;
plot(BW_list, log10(kappa_n), 's-', 'LineWidth', 1.2);
ylabel('log_{10}(条件数 \kappa(H_n))');
yyaxis right;
plot(BW_list, 100 * p_lin, 'o-', 'LineWidth', 1.2);
ylabel('线性分量功率占比 (%)');
xlabel('采集带宽 (MHz)');
set(gca, 'XDir', 'reverse');
grid on;
title('带宽降低：条件数升高、线性占比上升');
