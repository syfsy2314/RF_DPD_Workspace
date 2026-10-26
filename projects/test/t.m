%用于测试随着滤波器带宽的降低，条件数的变化（带宽扫描版）
%输出两条曲线：cond(H_x)（LS 实际承受的病态）与 cond(H_n)（列归一化，纯共线性）
%%
clear;
close all;
fs = 1000;
%% 导入数据（与带宽无关，只加载一次）
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
%% 扫描设置
BW_list = 100:-5:30;             % 采集带宽扫描（MHz），从大到小
nBW = numel(BW_list);
kappa_x = zeros(1, nBW);         % 原始基条件数（LS 问题口径）
kappa_n = zeros(1, nBW);         % 列归一化基条件数（纯共线性口径）
K = 7;
M = 5;
Dstop1 = 0.0001;
Dstop2 = 0.000001;
dens  = 20;
Dpass1 = (10^(0.5/20) - 1) / (10^(0.5/20) + 1);      % H1 通带波纹 0.5 dB
Dpass2 = 0.0057501127785;                            % H2 通带波纹（原脚本取值）
%% 带宽扫描
for iBW = 1:nBW
    BW = BW_list(iBW);
    % 第一个 FIR 滤波器（OMUX 等效，通带随 BW 缩放）
    [N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(fs/2), [1 0], [Dpass1, Dstop1]);
    H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
    % 第二个 FIR 滤波器（反馈覆盖）
    [N2, Fo, Ao, W] = firpmord([BW, BW+3]/(fs/2), [1 0], [Dpass2, Dstop2]);
    H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
    h2 = H2_FIR.Numerator;
    % 观测链路（对齐结果存独立变量，避免污染原始 x）
    y_f = H1_FIR(y0);
    [y_f, x_al, ~] = align_circular_signal(y_f, x);
    % 带限基函数及其条件数
    H_x = geth_even_BLMP(x_al, K, M, h2);
    Hn = H_x ./ vecnorm(H_x);                        % 列归一化
    s1 = svd(H_x, 'econ');
    s2 = svd(Hn,  'econ');
    kappa_x(iBW) = s1(1) / s1(end);
    kappa_n(iBW) = s2(1) / s2(end);
    fprintf('BW = %3d MHz: cond(H_x) = %.3g, cond(H_n) = %.3g\n', ...
        BW, kappa_x(iBW), kappa_n(iBW));
end
%% 绘图（纵轴为 log10(条件数)，线性坐标）
figure;
plot(BW_list, log10(kappa_x), 'o-', 'LineWidth', 1.2); hold on;
plot(BW_list, log10(kappa_n), 's-', 'LineWidth', 1.2);
xlabel('采集带宽 (MHz)');
ylabel('log_{10}(条件数)');
set(gca, 'XDir', 'reverse');         % 横轴从大到小显示（90 -> 50）
legend('\kappa(H_x) 原始基', '\kappa(H_n) 列归一化基', 'Location', 'best');
title('BLMP 基函数矩阵条件数随采集带宽的变化');
grid on;
