%% MG 幅度算子单次交替迭代测试
%% 前置处理
clear;
close all;
fs = 1000;
nfft = 1024;
window = hann(nfft);

%% 第一个 FIR 滤波器
BW = 70;
Fs = 1000;
Fpass = BW + 5;
Fstop = BW + 10;
Rp_dB = 3.1;
Dpass = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
Dstop = 0.0001;
dens = 20;
[N, Fo, Ao, W] = firpmord( ...
    [Fpass, Fstop] / (Fs/2), [1 0], [Dpass, Dstop]);
b = firpm(N, Fo, Ao, W, {dens});
H1_FIR = dsp.FIRFilter('Numerator', b);
h1 = H1_FIR.Numerator;

%% 第二个 FIR 滤波器
Fpass = BW;
Fstop = BW + 3;
Dpass = 0.0057501127785;
Dstop = 0.0001;
[N, Fo, Ao, W] = firpmord( ...
    [Fpass, Fstop] / (Fs/2), [1 0], [Dpass, Dstop]);
b = firpm(N, Fo, Ao, W, {dens});
H2_FIR = dsp.FIRFilter('Numerator', b);
h2 = H2_FIR.Numerator;

%% 初始化
load C:\Users\admin\Desktop\测试\第四次\x.mat;
load C:\Users\admin\Desktop\测试\第四次\y0.mat;
gain_rmspa = rms(y0) / rms(x);
K = 7;
M = 5;
s = 3;

%% 模拟反馈链路
y_f = H1_FIR(y0);
[y_f, x, ~] = align_circular_signal(y_f, x);
y_bf = H2_FIR(y_f);

%% 性能上限对照组
h_Ref = conv(h1, h2);
linear_reference_Ref = filter(h_Ref, 1, x);
y_bf_Ref = align_circular_signal(y_bf, linear_reference_Ref);
[y_Ref_raw, ~] = get_MG(K, M, s, x, h_Ref, y_bf_Ref);
y_Ref = align_circular_signal(y_Ref_raw, y0);
gain_rms_Ref = rms(y_Ref) / rms(x);

%% 普通覆盖对照组，相当于第 0 次迭代
[y_Initial, gain_rms_Initial] = get_MG(K, M, s, x, h2, y_bf);

%% 固定进行一次外层交替迭代
operator_R = 250;
operator_eta = 1e-4;
h_Iteration1 = estimate_amplitude_operator( ...
    y_Initial, y_bf, h2, operator_R, operator_eta);

linear_reference_Iteration1 = filter(h_Iteration1, 1, x);
y_bf_Iteration1 = align_circular_signal( ...
    y_bf, linear_reference_Iteration1);
[y_Iteration1, gain_rms_Iteration1] = get_MG( ...
    K, M, s, x, h_Iteration1, y_bf_Iteration1);

% 第 1 列为初始结果，第 2 列为一次迭代结果。
y_IterationResults = [y_Initial, y_Iteration1];
gain_rms_IterationResults = [gain_rms_Initial; gain_rms_Iteration1];

%% 展示结果
figure(1);
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');

[p1, f1] = pwelch(y0 / gain_rmspa, window, nfft/2, nfft, fs);
plot(f1 - fs/2, 10 * log10(fftshift(p1 * 100)), 'r-');
hold on;

[p4, f4] = pwelch(y_Ref / gain_rms_Ref, window, nfft/2, nfft, fs);
plot(f4 - fs/2, 10 * log10(fftshift(p4 * 100)), 'b-');

[p4, f4] = pwelch( ...
    y_Initial / gain_rms_Initial, window, nfft/2, nfft, fs);
plot(f4 - fs/2, 10 * log10(fftshift(p4 * 100)), 'g-');

[p4, f4] = pwelch( ...
    y_Iteration1 / gain_rms_Iteration1, window, nfft/2, nfft, fs);
plot(f4 - fs/2, 10 * log10(fftshift(p4 * 100)), 'c-');

legend('PA输出信号', 'Ref', 'Initial', 'Iteration-1');
grid on;

NMSE_Ref = get_NMSE(y_Ref, y0);
fprintf('参考NMSE为 %.2f dB\n',NMSE_Ref);
NMSE_IterationResults = [ ...
    get_NMSE(y_Initial, y0); ...
    get_NMSE(y_Iteration1, y0)];
NMSE_Initial = NMSE_IterationResults(1);
fprintf('首次估计NMSE为 %.2f dB\n',NMSE_Initial);
NMSE_Iteration1 = NMSE_IterationResults(2);
fprintf('最终迭代估计NMSE为 %.2f dB\n',NMSE_Iteration1);
