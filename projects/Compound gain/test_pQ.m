% 用复增益作为先验信息 验证平坦度不一样的时候的结果√
% 验证相位非线性情况下，能否实现
%% 
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第 一 个IIR滤波器 OMUX
BW = 60;
Fpass = BW+5;              % Passband Frequency
Fstop = BW+10;              % Stopband Frequency 
Apass = 0.01;    
Astop = 60;    
Fs    = 1000;  
h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);
H1_FIR = design(h, 'ellip', ...
    'MatchExactly', 'both', ...
    'SystemObject', true);
h1 = H1_FIR.Numerator;
%% 第一个FIR滤波器
% BW = 70;
% Fs = 1000;  % Sampling Frequency
% Fpass = BW+5;              % Passband Frequency
% Fstop = BW+10;              % Stopband Frequency
% % Dpass = 0.57501127785;  % Passband Ripple
% Rp_dB = 3.1;        % 通带波纹，单位 dB 大概是0.17 %% 3dB有bug 3.1都行
% Dpass = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
% Dstop = 0.0001;          % Stopband Attenuation
% dens  = 20;              % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% b  = firpm(N, Fo, Ao, W, {dens});
% H1_FIR = dsp.FIRFilter( ...
%     'Numerator', b);
% h1 = H1_FIR.Numerator;
%% 第二个FIR滤波器
%UNTITLED 返回离散时间滤波器对象。
Fs = 1000; 
Fpass = BW;             
Fstop = BW+3;                
Dpass = 0.0057501127785;  
Dstop = 0.0001;          
dens  = 20;             
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H2_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h2 =  H2_FIR.Numerator;
%%
load C:\Users\admin\Desktop\测试\第四次\x.mat;
load C:\Users\admin\Desktop\测试\第四次\y0.mat; 
gain_rmspa = rms(y0)/rms(x);
%% 先进行 小算子的估计
x_f = H1_FIR(x);
[x_f,x,~] = align_circular_signal(x_f,x);
g = x./x_f;
%% OMUX模拟滤波
y_f = H1_FIR(y0);
[y_f,x,~] = align_circular_signal(y_f,x);
% FIR滤波
y_bf = H2_FIR(y_f .* g);
reset(H2_FIR);
h = h2;
%% 初始化
K = 7;
M = 5;
L = K*(M+1);
H_x = geth_even_BLMP(x,K,M,h);%中频0 带限基函数
X = geth_even_1d(x,K,M);
%% 初始估计
S = 0.3;
lamda = 1;
T_delta = 1e-4;
kappa_G_max = 1e4;
[y_Initial,~] = get_ALVS_LD_PLS(K,M,X,x,h,H2_FIR(y_f),S,lamda,T_delta,kappa_G_max);
% [y_Initial,~] = get_OPLS(K,M,X,x,h,H2_FIR(y_f),S,lamda);
%% 信号恢复
s = 2;
% MG
[y_MG,gain_rmsMG] = get_MG(K,M,s,x,h,y_bf);
%%
[y_LD_PLS,gain_rmsOPLS] = get_ALVS_LD_PLS(K,M,X,x,h,y_bf,S,lamda,T_delta,kappa_G_max);
% [y_LD_PLS,~] = get_OPLS(K,M,X,x,h,y_bf,S,lamda);
%% 
fprintf('初始估计NMSE为 %.2f dB\n',get_NMSE(y_Initial,y0));
fprintf('MG估计NMSE为 %.2f dB\n',get_NMSE(y_MG,y0));
fprintf('LD_PLS估计NMSE为 %.2f dB\n',get_NMSE( y_LD_PLS,y0));


