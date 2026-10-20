%% 用线性响应估计复增益
%% 前置处理
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第 一 个IIR滤波器 OMUX
% Fpass = 75;    % Passband Frequency
% Fstop = 85;    % Stopband Frequency
% Apass = 0.01;     % Passband Ripple (dB)
% Astop = 60;    % Stopband Attenuation (dB)
% Fs    = 1000;  % Sampling Frequency
% h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);
% Hd_OMUX = design(h, 'ellip', ...
%     'MatchExactly', 'both', ...
%     'SystemObject', true);
%% 第一个FIR滤波器
%UNTITLED 返回离散时间滤波器对象。
BW = 70;
Fs = 1000;  % Sampling Frequency
Fpass = BW+5;              % Passband Frequency
Fstop = BW+10;              % Stopband Frequency
% Dpass = 0.57501127785;  % Passband Ripple
Rp_dB = 3.1;        % 通带波纹，单位 dB 大概是0.17 %% 3dB有bug 3.1都行
Dpass = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H1_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h1 = H1_FIR.Numerator;
% [EOF]
Hd_OMUX = H1_FIR;
%% 第二个FIR滤波器
%UNTITLED 返回离散时间滤波器对象。
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
%%
load C:\Users\admin\Desktop\测试\第四次\x.mat;
load C:\Users\admin\Desktop\测试\第四次\y0.mat;
gain_rmspa = rms(y0)/rms(x);
N = length(x);
%% ======================OMUX模拟滤波========================
y_f = Hd_OMUX(y0);
% 时间对齐
[y_f,x,~] = align_circular_signal(y_f,x);
h = h2;
y_bf = H2_FIR(y_f);
%% ======================初始化=============================
K = 7;
M = 5;
L = K*(M+1);
H_x = geth_even_BLMP(x,K,M,h);%中频0 带限基函数
X = geth_even_1d(x,K,M);
%% =======================对照组=============================
s = 3;
% MG
[y_MG,gain_rmsMG] = get_MG(K,M,s,x,h,y_bf);
NMSE_MG = get_NMSE(y_MG, y0);
%% =====================分离并计算线性响应=================
Len = 1;
x_lin = H_x(:,Len);
X_nl_raw = H_x(:,(Len+1):end);
w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);
%计算线性响应
y1 = x_lin * w_lin;
gain_rmsL = rms(y1)/rms(x);
% 
figure(2);
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
[p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
hold on;
[p1, f1] = pwelch(x, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'b-')%红色
hold on;
[p4, f4] = pwelch(y1/gain_rmsL, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'g-')
hold on;
legend('PA输出信号','输入信号','线性');
grid on;
%% =====================估计复增益==========================================
% x_f = H1_FIR(x);
% [x_f,x,~] = align_circular_signal(x_f,x);
% g = x./x_f;
%% =====================线性响应估计复增益================================
[y1,x,~] = align_circular_signal(y1,x);
g = x./y1;
%% =====================复增益FIR滤波============================================
y_bf = H2_FIR(y_f .* g);
%% =====================信号恢复 第二轮迭代===========================================
s = 3;
% MG
[y_MG_propose,gain_rmsMG_propose] = get_MG(K,M,s,x,h,y_bf);
NMSE_propose = get_NMSE(y_MG_propose, y0);
%% =====================展示结果=======================================
figure(1);
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
[p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
hold on;
[p4, f4] = pwelch(y_MG_propose/gain_rmsMG_propose, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'b-')
hold on;
[p4, f4] = pwelch(y_MG/gain_rmsMG, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'g-')
hold on;
legend('PA输出信号','MG-propose','MG');
grid on;
%%
S = 0.1;
lamda = 0.01;
[y_OPLS,gain_rmsOPLS] = get_OPLS(K,M,X,x,h,y_bf,S,lamda);
fprintf('第一次计算的NMSE = %.2f dB\n',get_NMSE(y_OPLS,y0));
fprintf('===========================================================\n');