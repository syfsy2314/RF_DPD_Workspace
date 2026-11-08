%% 忽略第一个FIR 平坦情况 不迭代
%%
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第一个FIR滤波器
BW = 50;
%UNTITLED 返回离散时间滤波器对象。
Fs = 1000;  % Sampling Frequency
Fpass = BW+5;              % Passband Frequency
Fstop = BW+10;              % Stopband Frequency
% Dpass = 0.57501127785;  % Passband Ripple
Rp_dB = 0.5;        % 通带波纹，单位 3dB 大概是0.17
Dpass = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H1_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h1 = H1_FIR.Numerator;
% [EOF]
%% 第二个FIR滤波器
%UNTITLED 返回离散时间滤波器对象。
Fs = 1000; 
Fpass = BW;             
Fstop = BW+3;              
Dpass = 0.0057501127785;
Dstop = 1e-6;            % 数字域覆盖滤波器，深阻带排除带外泄漏
dens  = 20;             
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H2_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h2 =  H2_FIR.Numerator;
%% 导入数据
load E:\测试\2026.6.10\test1\x_v0.6_v1.mat;
load E:\测试\2026.6.10\test1\y0_v0.6_v1.mat;
gain_rmspa = rms(y0)/rms(x);
N = length(x);
%% 模拟带限滤波器
y_f = H1_FIR(y0);
[y_f,x,~] = align_circular_signal(y_f,x);
% FIR滤波
y_bf = H2_FIR(y_f);
% h = conv(h1,h2);
h = h2;
%% 初始化
%================================
% MP
K = 7;
M = 5;
L = K*(M+1);
H_x = geth_even_BLMP(x,K,M,h);%中频0 带限基函数
X = geth_even_1d(x,K,M);
%=================================
% a = inv(ctranspose(H_x) * H_x) * ctranspose(H_x) * y_bf;
% y_BLMP = X * a;
% NESE_BLMP = get_NMSE(y_BLMP,y0);
%=================================
%% MG
s = 2;
% MG
[y_MG,gain_rmsMG] = get_MG(K,M,s,x,h,y_bf);
fprintf('MG计算完成\n===================================================\n');
%% MR
Kb = 2;
nb = 0.2;
a1 = 0.0001;
Lb = get_Lb(H_x,Kb,nb);
%
Kg = 20;
ng = 0.2;
a2 = 0.00001;
g = 1000;
G = H_x(1:g,:);
Lg = get_Lg(G,Kg,ng);
%
a3 = 0.1;
%
pa2= inv(ctranspose(H_x) * H_x + a1 * Lb + a2 * G' * Lg * G) * ctranspose(H_x) * y_bf;
y_MR = X * pa2;
fprintf('MR计算完成\n===================================================\n');
gain_rmsMR = rms(y_MR)/rms(x);
%% SBL
[y_SBL,gain_rmsSBL] = get_SBL(K,M,X,H_x,x,y_bf,3);
fprintf('SBL计算完成\n===================================================\n');
%%
S = 0.3;
lamda = 10;
[y_OPLS,gain_rmsOPLS] = get_OPLS(K,M,X,x,h,y_bf,S,lamda);
fprintf('LD-PLS计算完成\n===================================================\n');
% %% 第二次
% hh = ls_filter_est(y_OPLS, y_bf, 2000);
% [y_OPLS2,gain_rmsOPLS2] = get_OPLS(K,M,X,x,hh,y_bf,S,lamda);
% %% 第三次
% hhh = ls_filter_est(y_OPLS2, y_bf, 2000);
% [y_OPLS3,gain_rmsOPLS3] = get_OPLS(K,M,X,x,hhh,y_bf,S,lamda);
%% 绘图
figure(1);
[p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
hold on;
[p4, f4] = pwelch(y_MG/gain_rmsMG, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'g-')
hold on;
[p4, f4] = pwelch(y_OPLS/gain_rmsOPLS, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'm-')
hold on;
[p4, f4] = pwelch(y_MR/gain_rmsMR, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'k-')
hold on;
[p4, f4] = pwelch(y_SBL/gain_rmsSBL, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'b-')
hold on;
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
legend('PA输出信号','MG','OPLS','MR','SBL');
grid on;
%% 指标
NESE_MG = get_NMSE(y_MG,y0);
NESE_LD_PLS = get_NMSE(y_OPLS,y0);
NESE_MR = get_NMSE(y_MR,y0);
NESE_SBL = get_NMSE(y_SBL,y0);
% NESE_OPLS2 = get_nmse(y_OPLS2,y0);
% NESE_OPLS3 = get_nmse(y_OPLS3,y0);

