%% 估计低维幅度算子
%% 前置处理
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第一个FIR滤波器
BW = 70;
Fs = 1000;  % Sampling Frequency
Fpass = BW+5;              % Passband Frequency
Fstop = BW+10;              % Stopband Frequency
% Dpass = 0.57501127785;  % Passband Ripple
Rp_dB = 5.8;        % 通带波纹，单位 dB 大概是0.17 %% 3dB有bug 3.1都行
Dpass = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H1_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h1 = H1_FIR.Numerator;
%% 第二个FIR滤波器
Fs = 1000;  % Sampling Frequency
Fpass = BW;              % Passband Frequency
Fstop = BW+3;              % Stopband Frequency
Dpass = 0.0057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H2_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h2 =  H2_FIR.Numerator;
%% 初始化
load C:\Users\admin\Desktop\测试\第四次\x.mat;
load C:\Users\admin\Desktop\测试\第四次\y0.mat;
gain_rmspa = rms(y0)/rms(x);
N = length(x);
K = 7; M = 5;
S = 0.3;
lamda = 1;
X = geth_even_1d(x, K, M);
%% ======================模拟滤波========================
y_f = H1_FIR(y0);
% 时间对齐
[y_f,x,~] = align_circular_signal(y_f,x);
h = h2;
y_bf = H2_FIR(y_f);
%% =======================性能上限对照组==========================
h_Ref = conv(h1, h2);
% h1 和 h2 均按因果 FIR 实现，会在模型侧引入群时延。
% 将 y_bf 的副本统一对齐到复合滤波器的线性参考，避免反馈样本与
% 带限 BLMP 回归矩阵发生时间错位；原始 y_bf 保持不变。
linear_reference_Ref = filter(h_Ref, 1, x);
[y_bf_Ref, ~, delay_feedback_Ref] = ...
    align_circular_signal(y_bf, linear_reference_Ref);
% 使用完全匹配的反馈算子估计 PA 参数并恢复全带输出。
[y_Ref_raw, ~] = get_OPLS(K, M, X, x, h_Ref, y_bf_Ref, 0.3, lamda);
% 与真实 y0 的对齐只用于离线性能评价，不参与上限组的参数估计。
[y_Ref, ~, delay_output_Ref] = align_circular_signal(y_Ref_raw, y0);
gain_rms_Ref = rms(y_Ref) / rms(x);
%% =======================普通覆盖对照组 相当于迭代 0 次==========================
[y_Initial, gain_rms_Initial] = get_OPLS(K, M, X, x, h2, y_bf, S, lamda);
%% ========================交替迭代===============================
% iteration_count 为外层交替迭代次数，可根据实验需要修改。
% y_IterationResults 的第 1 列保存初始恢复（第 0 次迭代），
% 第 i+1 列保存第 i 次外层迭代结果，确保历史输出不会被覆盖。
iteration_count = 1;
operator_R = 250;
operator_eta = 1e-4;
[y_IterationResults, gain_rms_IterationResults] = ...
    run_amplitude_alternating_iterations( ...
    iteration_count, y_Initial, gain_rms_Initial, y_bf, h2, ...
    x, K, M, X, S, lamda, operator_R, operator_eta);
% 最后一列为当前设置下的最终恢复结果。
y_IterationFinal = y_IterationResults(:, end);
gain_rms_IterationFinal = gain_rms_IterationResults(end);
%% =========================展示结果========================================
% figure(1);
% xlabel('频率偏移 (MHz)');
% ylabel('功率谱密度 (dBm/Hz)');
% [p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
% plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
% hold on;
% [p4, f4] = pwelch(y_Ref/gain_rms_Ref, window, nfft/2, nfft ,fs);
% plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'b-')
% hold on;
% [p4, f4] = pwelch(y_Initial/gain_rms_Initial, window, nfft/2, nfft ,fs);
% plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'g-')
% hold on;
% [p4, f4] = pwelch(y_IterationFinal / gain_rms_IterationFinal, ...
%     window, nfft/2, nfft, fs);
% plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'c-')
% legend('PA输出信号', 'Ref', 'Initial', ...
%     sprintf('Iteration-%d', iteration_count));
% grid on;
% NMSE_IterationResults 与输出矩阵列一一对应：第 1 个元素为初始结果。
NMSE_IterationResults = zeros(iteration_count + 1, 1);
for result_index = 1:(iteration_count + 1)
    NMSE_IterationResults(result_index) = get_NMSE( ...
        y_IterationResults(:, result_index), y0);
end
NMSE_Initial = NMSE_IterationResults(1);
NMSE_IterationFinal = NMSE_IterationResults(end);
NMSE_Ref = get_NMSE(y_Ref, y0);
fprintf('=====================================================\n');
fprintf('参考NMSE为 %.2f dB\n',NMSE_Ref);
fprintf('首次估计NMSE为 %.2f dB\n',NMSE_Initial);
fprintf('最终迭代估计NMSE为 %.2f dB\n',NMSE_IterationFinal);
