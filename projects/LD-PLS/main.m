%% 忽略第一个FIR 平坦情况 不迭代
%%
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第一个FIR滤波器
BW = 91;
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
projectDir = fileparts(mfilename('fullpath'));
datasetsDir = fullfile(fileparts(fileparts(projectDir)), 'datasets');
test2Dir = fullfile(datasetsDir, '测试', '2026.11.24', 'test2');
dataY = load(fullfile(test2Dir, 'y0.mat'), 'y0');
y0 = dataY.y0(:);
[~, ~, x] = import_iq_binary(fullfile(datasetsDir, 'data_38400_v0.6_test1.dat'), numel(y0));
assert(numel(x) == numel(y0), 'test2 输入与参考输出的样点数必须一致。');
x = x / max(abs(x));
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
s = 3;
% MG
[y_MG,gain_rmsMG] = get_MG(K,M,s,x,h,y_bf);
fprintf('MG                         %8.3f dB\n', get_NMSE(y_MG, y0));
fprintf('MG计算完成\n===================================================\n');
%% MR
Kb = 2;
nb = 0.2;
a1 = 1e-4;
Lb = get_Lb(H_x,Kb,nb);
%
Kg = 20;
ng = 0.2;
a2 = 0.000000001;
g = 1000;
G = H_x(1:g,:);
Lg = get_Lg(G,Kg,ng);
%
a3 = 0.1;
%
pa2= inv(ctranspose(H_x) * H_x + a1 * Lb + a2 * G' * Lg * G) * ctranspose(H_x) * y_bf;
y_MR = X * pa2;
fprintf('MR                         %8.3f dB\n', get_NMSE(y_MR, y0));
fprintf('MR计算完成\n===================================================\n');
gain_rmsMR = rms(y_MR)/rms(x);
%% SBL
[y_SBL,gain_rmsSBL] = get_SBL(K,M,X,H_x,x,y_bf,3);
fprintf('SBL计算完成\n===================================================\n');
%%
lamda = 10;
T_delta = 0;
kappa_G_max = 10^4.5;
[y_LD_PLS,gain_rmsLD_PLS,O] = get_adaptive_LD_PLS( ...
    K,M,X,x,h,y_bf,lamda,T_delta,kappa_G_max);
[hFullNoPenalty, O_no_penalty] = get_energy_regularized_LD_PLS( ...
    X, H_x, y_bf, lamda, 0, T_delta, kappa_G_max, size(X, 2) - 1);
y_LD_PLS_no_penalty = X * hFullNoPenalty;
gain_rmsLD_PLS_no_penalty = rms(y_LD_PLS_no_penalty) / rms(x);
fprintf('全带宽能量惩罚 LD-PLS 计算完成，O=%d\n===================================================\n', O);
% 固定潜变量数时可设 S=0.3 并调用 get_LD_PLS(K,M,X,x,h,y_bf,S,lamda)。
% %% 第二次
% hh = ls_filter_est(y_LD_PLS, y_bf, 2000);
% [y_LD_PLS2,gain_rmsLD_PLS2] = get_LD_PLS(K,M,X,x,hh,y_bf,S,lamda);
% %% 第三次
% hhh = ls_filter_est(y_LD_PLS2, y_bf, 2000);
% [y_LD_PLS3,gain_rmsLD_PLS3] = get_LD_PLS(K,M,X,x,hhh,y_bf,S,lamda);
%% 全带 NMSE
fprintf('\n全带 NMSE（观测半带宽 %g MHz）\n', BW);

fprintf('有全带惩罚 LD-PLS          %8.3f dB\n', get_NMSE(y_LD_PLS, y0));
fprintf('无全带惩罚 LD-PLS          %8.3f dB\n', get_NMSE(y_LD_PLS_no_penalty, y0));

fprintf('SBL                        %8.3f dB\n', get_NMSE(y_SBL, y0));
%% 绘图
figure(1);
[p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
hold on;
[p4, f4] = pwelch(y_MG/gain_rmsMG, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'g-')
hold on;
[p4, f4] = pwelch(y_LD_PLS/gain_rmsLD_PLS, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'm-')
hold on;
[p4, f4] = pwelch(y_LD_PLS_no_penalty / gain_rmsLD_PLS_no_penalty, window, nfft/2, nfft, fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'm--')
hold on;
[p4, f4] = pwelch(y_MR/gain_rmsMR, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'k-')
hold on;
[p4, f4] = pwelch(y_SBL/gain_rmsSBL, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'b-')
hold on;
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
legend('PA输出信号','MG','LD-PLS (\mu=1)','LD-PLS (\mu=0)','MR','SBL');
grid on;

