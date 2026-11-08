%% 幅度自校准：通带纹波扫描（0.1 -> 2 dB，步进 0.2 dB）
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
Rp_list = 0.1:0.2:2;             % H1 通带纹波扫描（dB）
nRp = numel(Rp_list);
NMSE_Ref  = zeros(nRp, 1);
NMSE_Init = zeros(nRp, 1);
NMSE_Iter = zeros(nRp, 1);
%% 固定参数
BW = 80;
Fs = 1000;
Dpass2 = 0.0057501127785;        % H2 通带波纹（数字域，固定）
Dstop2 = 1e-4;                 % H2 阻带（数字域，固定）
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
%% 纹波扫描（H1 固定阶数：平坦低通 × 独立通带余弦波纹）
% 平坦低通：一次设计、阶数固定；通带纹波 0.01 dB，远小于最小扫描 Rp（0.1 dB）
Dpass_flat = (10^(0.01/20) - 1) / (10^(0.01/20) + 1);
[N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(Fs/2), [1 0], [Dpass_flat, 1e-4]);
h_flat = firpm(N1, Fo, Ao, W, {dens});
nRipple = 5;                                    % 通带内波纹周期数（固定）
tau_ripple = round(nRipple * Fs / (BW+5));      % 波纹周期对应的样本间隔
Rp_actual = zeros(nRp, 1);      % 实测通带纹波（dB）
for iRp = 1:nRp
    Rp_dB = Rp_list(iRp);
    % ---- H1 = h_flat ⊛ h_ripple（频域幅度相乘 ⟷ 时域卷积） ----
    % 余弦幅度波纹 1+ε·cos(2πf·τ/Fs) 的精确时域表示：3 抽头零相位 FIR [ε/2, 1, ε/2]，
    % 无逼近误差；h_ripple 群延迟为常数 τ，与 h_flat 一起被对齐吸收
    epsilon = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);   % 波纹线性幅度，峰峰 = Rp_dB
    h_ripple = [epsilon/2, zeros(1, tau_ripple-1), 1, zeros(1, tau_ripple-1), epsilon/2];
    h1 = conv(h_flat, h_ripple);
    H1_FIR = dsp.FIRFilter('Numerator', h1);
    % 实测通带纹波
    [Hmag, fgrid] = freqz(h1, 1, 4096, Fs);
    passband = fgrid <= (BW + 5);
    Rp_actual(iRp) = 20 * log10(max(abs(Hmag(passband))) / min(abs(Hmag(passband))));
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
    fprintf('纹波 标称%4.2f/实际%5.3f dB: Ref %7.2f | Initial %7.2f | 自校准 %7.2f dB\n', ...
        Rp_dB, Rp_actual(iRp), NMSE_Ref(iRp), NMSE_Init(iRp), NMSE_Iter(iRp));
end
%% 汇总表
fprintf('\n================ 纹波扫描汇总 ================\n');
fprintf('%8s %8s %12s %14s %14s\n', '标称(dB)', '实际(dB)', '参考NMSE', '首次估计NMSE', '自校准NMSE');
for iRp = 1:nRp
    fprintf('%8.2f %8.3f %12.2f %14.2f %14.2f\n', ...
        Rp_list(iRp), Rp_actual(iRp), NMSE_Ref(iRp), NMSE_Init(iRp), NMSE_Iter(iRp));
end
%% 曲线（横轴用标称纹波，阶数固定后连续无跳变）
figure;
plot(Rp_list, NMSE_Ref, 'b-o', 'LineWidth', 1.2); hold on;
plot(Rp_list, NMSE_Init, 'g-s', 'LineWidth', 1.2);
plot(Rp_list, NMSE_Iter, 'm-^', 'LineWidth', 1.2);
xlabel('H1 通带纹波 (dB)');
ylabel('NMSE (dB)');
legend('Ref 完美已知滤波器', 'Initial 传统覆盖', '一次自校准', 'Location', 'best');
title('幅度自校准性能随反馈链路纹波的变化');
grid on;
