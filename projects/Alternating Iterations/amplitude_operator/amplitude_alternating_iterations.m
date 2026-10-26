%% 幅度自校准：通带纹波扫描（0.01 -> 6 dB，步进 0.5 dB）
% 每个纹波点依次运行三组对照：
%   Ref      —— 完美已知复合滤波器 conv(h1,h2) 的性能上限
%   Initial  —— 传统覆盖（仅 h2，无视幅度失配），即未校准基线
%   Iteration—— 一次结构化幅度自校准后的恢复
% 输出：逐点打印 + 汇总表 + NMSE-纹波曲线
%%
clear;
close all;
fs = 1000;
nfft = 1024;
window = hann(nfft);
%% 扫描设置
Rp_list = 0.02:0.3:6;            % H1 通带纹波扫描（dB）
nRp = numel(Rp_list);
NMSE_Ref  = zeros(nRp, 1);
NMSE_Init = zeros(nRp, 1);
NMSE_Iter = zeros(nRp, 1);
%% 固定参数
BW = 80;
Fs = 1000;
Dpass2 = 0.0057501127785;        % H2 通带波纹（数字域，固定）
Dstop2 = 1e-4;                 % H2 阻带（数字域，固定）
Dstop1 = 0.0001;                 % H1 阻带（模拟域，固定）
dens = 20;
K = 7; M = 5;
S = 0.3;
lamda = 1;
iteration_count = 1;
operator_R = 250;
operator_eta = 1e-4;
%% 数据（与纹波无关，只加载一次）
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
gain_rmspa = rms(y0)/rms(x);
%% 纹波扫描
for iRp = 1:nRp
    Rp_dB = Rp_list(iRp);
    % ---- H1：模拟滤波器，通带纹波为扫描变量 ----
    Dpass1 = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
    [N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(Fs/2), [1 0], [Dpass1, Dstop1]);
    H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
    % ---- H2：数字覆盖滤波器（固定） ----
    [N2, Fo, Ao, W] = firpmord([BW, BW+3]/(Fs/2), [1 0], [Dpass2, Dstop2]);
    H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
    h2 = H2_FIR.Numerator;
    % ---- 观测链路（对齐结果存独立变量，避免污染原始 x） ----
    y_f = H1_FIR(y0);
    [y_f, x_al, ~] = align_circular_signal(y_f, x);
    y_bf = H2_FIR(y_f);
    X = geth_even_1d(x_al, K, M);
    % ---- Ref：完美已知复合滤波器（性能上限） ----
    h1 = H1_FIR.Numerator;
    h_Ref = conv(h1, h2);
    linear_reference_Ref = filter(h_Ref, 1, x_al);
    [y_bf_Ref, ~, ~] = align_circular_signal(y_bf, linear_reference_Ref);
    [y_Ref_raw, ~] = get_OPLS(K, M, X, x_al, h_Ref, y_bf_Ref, S, lamda);
    [y_Ref, ~, ~] = align_circular_signal(y_Ref_raw, y0);
    NMSE_Ref(iRp) = get_NMSE(y_Ref, y0);
    % ---- Initial：传统覆盖（迭代 0 次） ----
    [y_Initial, gain_rms_Initial] = get_OPLS(K, M, X, x_al, h2, y_bf, S, lamda);
    % ---- 一次结构化幅度自校准 ----
    [y_IterationResults, ~] = run_amplitude_alternating_iterations( ...
        iteration_count, y_Initial, gain_rms_Initial, y_bf, h2, ...
        x_al, K, M, X, S, lamda, operator_R, operator_eta);
    NMSE_Init(iRp) = get_NMSE(y_Initial, y0);
    NMSE_Iter(iRp) = get_NMSE(y_IterationResults(:, end), y0);
    fprintf('纹波 %4.2f dB: Ref %7.2f | Initial %7.2f | 自校准 %7.2f dB\n', ...
        Rp_dB, NMSE_Ref(iRp), NMSE_Init(iRp), NMSE_Iter(iRp));
end
%% 汇总表
fprintf('\n================ 纹波扫描汇总 ================\n');
fprintf('%8s %12s %14s %14s\n', '纹波(dB)', '参考NMSE', '首次估计NMSE', '自校准NMSE');
for iRp = 1:nRp
    fprintf('%8.2f %12.2f %14.2f %14.2f\n', ...
        Rp_list(iRp), NMSE_Ref(iRp), NMSE_Init(iRp), NMSE_Iter(iRp));
end
%% 曲线
figure;
plot(Rp_list, NMSE_Ref, 'b-o', 'LineWidth', 1.2); hold on;
plot(Rp_list, NMSE_Init, 'g-s', 'LineWidth', 1.2);
plot(Rp_list, NMSE_Iter, 'm-^', 'LineWidth', 1.2);
xlabel('H1 通带纹波 (dB)');
ylabel('NMSE (dB)');
legend('Ref 完美已知滤波器', 'Initial 传统覆盖', '一次自校准', 'Location', 'best');
title('幅度自校准性能随反馈链路纹波的变化');
grid on;
