%% 相位校准（Phase Calibration）项目入口
% H1 用椭圆 IIR 滤波器模拟实际模拟滤波器的非线性相位（群延迟变化）。
% 流程：
%   1) 简单覆盖（只用 h2 FIR）—— 基线，忽略相位，性能差
%   2) 相位校准（非对称算子盲估幅度+相位）—— 从信号盲估，不需要 H1 信息
clear; close all;
%% 参数
BW = 80;
Fs = 1000;
Rp_dB = 0.5;                 % H1 通带纹波（椭圆 IIR）
Rs_dB = 80;                  % H1 阻带衰减（dB）
Dpass2 = 0.0057501127785;    % H2 通带
Dstop2 = 1e-4;               % H2 阻带
dens = 20;
K = 7; M = 5;
S = 0.3;
lamda = 1;
N_fir = 3000;                % IIR 冲激响应截断长度
operator_R = 1000;           % 非对称算子单边阶数（完整 2R+1）
operator_eta = 1e-4;
%% 数据
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
%% H1：椭圆 IIR 滤波器（含非线性相位）
[N1, Wn] = ellipord((BW+5)/(Fs/2), (BW+10)/(Fs/2), Rp_dB, Rs_dB);
[b1, a1] = ellip(N1, Rp_dB, Rs_dB, Wn, 'low');
fprintf('H1(椭圆 IIR) 阶数 = %d\n', N1);
h1_fir = impz(b1, a1, N_fir);   % 冲激响应截断成 FIR，保留幅度+相位
%% H2：FIR（数字覆盖）
[N2, Fo, Ao, W] = firpmord([BW, BW+3]/(Fs/2), [1 0], [Dpass2, Dstop2]);
h2 = firpm(N2, Fo, Ao, W, {dens});
%% 观测链路
y_f = filter(h1_fir, 1, y0);     % H1（椭圆响应，FIR 实现）
[y_f, x_al, ~] = align_circular_signal(y_f, x);
y_bf = filter(h2, 1, y_f);       % H2（FIR 覆盖）
X = geth_even_1d(x_al, K, M);
%% 基线：简单覆盖（只用 h2，忽略相位）
[y_cover, ~] = get_OPLS(K, M, X, x_al, h2, y_bf, S, lamda);
NESE_cover = get_NMSE(y_cover, y0);
%% 相位校准：非对称算子盲估 + 交替迭代（让 yCoarse 逐步逼近 y0）
iteration_count = 5;
yIter = zeros(length(y_cover), iteration_count + 1);
yIter(:, 1) = y_cover;
fprintf('迭代 0: NMSE = %.2f dB\n', NESE_cover);
for i = 1:iteration_count
    [hModel, ~, ~, fitDb] = estimate_complex_operator( ...
        yIter(:, i), y_bf, h2, operator_R, operator_eta);
    linearReference = filter(hModel, 1, x_al);
    [y_bf_aligned, ~, ~] = align_circular_signal(y_bf, linearReference);
    [yIter(:, i+1), ~] = get_OPLS(K, M, X, x_al, hModel, y_bf_aligned, S, lamda);
    fprintf('迭代 %d: 拟合 %.2f dB, 建模 NMSE = %.2f dB\n', ...
        i, fitDb, get_NMSE(yIter(:, i+1), y0));
end
NESE_cal = get_NMSE(yIter(:, end), y0);
%% 结果
fprintf('================================================\n');
fprintf('简单覆盖 NMSE = %.2f dB\n', NESE_cover);
fprintf('相位校准 NMSE = %.2f dB（迭代 %d 次）\n', NESE_cal, iteration_count);
fprintf('改善           = %.2f dB\n', NESE_cal - NESE_cover);
fprintf('================================================\n');
