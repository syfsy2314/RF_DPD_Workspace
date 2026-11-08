%%验证DOMP算法的性能
%%
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第一个FIR滤波器
BW = 100;
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
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
gain_rmspa = rms(y0)/rms(x);
N = length(x);
%% OMUX模拟滤波
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
%%
Rmat_label = cell(1, K*(M+1));   % 42 个空 cell，打印时核名留空，其余正常
[h_domp, s, nopt, h_full] = omp_domp(H_x, y_bf, Rmat_label, 'DOMP');
% BIC 在大样本+带限病态场景下惩罚不足，会选入"暗方向"核导致带外重建崩坏，
% 改用全带重建 NMSE 扫描定阶
nmse_rec = zeros(1, K*(M+1));
for n = 1:K*(M+1)
    nmse_rec(n) = get_NMSE(X*h_full(:,n), y0);
end
[~, n_best] = min(nmse_rec);
y_domp = X * h_full(:, n_best);      % 全带基重建，与 MR/SBL/OPLS 同口径
NESE_DOMP = get_NMSE(y_domp, y0);    % 与 main.m 指标段并列
fprintf('DOMP: BIC定阶 nopt=%d -> 重建 %.2f dB\n', nopt, get_NMSE(X*h_domp, y0));
fprintf('DOMP: 扫描定阶 n=%d -> 重建 %.2f dB\n', n_best, NESE_DOMP);
% figure;
% plot(nmse_rec); hold on;
% plot(n_best, nmse_rec(n_best), 'ro');
% xlabel('核个数'); ylabel('重建 NMSE (dB)'); title('DOMP 定阶扫描');
% grid on;
%% OPLS
% S = 0.3;
% lamda = 10;
% [y_LD_PLS,gain_rmsOPLS] = get_OPLS(K,M,X,x,h,y_bf,S,lamda);
% NESE_LD_PLS = get_NMSE(y_LD_PLS, y0);    
% fprintf('NMSE_LD_PLS = %.2f dB\n', NESE_LD_PLS);
% LD-PLS 扫描定阶(与 DOMP 扫描同口径): 扫 O 与 lamda, 判据为全带重建 NMSE
% 增量实现: 潜变量一次提取, 逐 O 复用(见 get_OPLS_full)
% 经验上限: O 不超过 0.5*T, 扫描范围取一半
T_all = K*(M+1);
O_scan = round(0.5*T_all);
lam_list = [0.1 1 10 100 1000];
nmse_opls = zeros(O_scan, numel(lam_list));    % O = 1..O_scan
for il = 1:numel(lam_list)
    Y_full = get_OPLS_full(K,M,X,x,h,y_bf,lam_list(il),O_scan);
    for O = 1:O_scan
        nmse_opls(O,il) = get_NMSE(Y_full(:,O), y0);
    end
    % if lam_list(il) == lamda
    %     % 一致性检查: 增量版与 get_OPLS(固定规则) 在同 O 处应严格一致
    %     fprintf('一致性检查(O=%d, 与get_OPLS最大偏差): %.2e\n', ...
    %         round(S*T_all), max(abs(Y_full(:,round(S*T_all)) - y_LD_PLS)));
    % end
end
[nmse_gmin, idx] = min(nmse_opls(:));
[O_best, il_best] = ind2sub(size(nmse_opls), idx);
% fprintf('OPLS: 固定规则 O=%d -> 重建 %.2f dB\n', round(S*T_all), NESE_LD_PLS);
% fprintf('OPLS: 扫描最优 O=%d, lamda=%g -> 重建 %.2f dB\n', ...
%     O_best, lam_list(il_best), nmse_gmin);
fprintf('NMSE_LD_PLS = %.2f dB\n', nmse_gmin);
% figure;
% plot(1:O_scan, nmse_opls, 'LineWidth', 1); hold on;
% plot(O_best, nmse_gmin, 'ro', 'MarkerFaceColor', 'r');
% xlabel('潜变量个数 O'); ylabel('重建 NMSE (dB)'); title('LD-PLS 定阶扫描');
% legend(arrayfun(@(l) sprintf('lamda=%g', l), lam_list, 'UniformOutput', false), ...
%     'Location', 'best');
% grid on;
%% MG
s = 3;
% MG
[y_MG,gain_rmsMG] = get_MG(K,M,s,x,h,y_bf);
NESE_MG = get_NMSE(y_MG, y0);    
fprintf('NMSE_MG = %.2f dB\n', NESE_MG);
%% MR
Kb = 2;
nb = 0.2;
a1 = 1e-4;
Lb = get_Lb(H_x,Kb,nb);
%
Kg = 20;
ng = 0.2;
a2 = 1e-7;
g = 1000;
G = H_x(1:g,:);
Lg = get_Lg(G,Kg,ng);
%
a3 = 0.1;
%
pa2= inv(ctranspose(H_x) * H_x + a1 * Lb + a2 * G' * Lg * G) * ctranspose(H_x) * y_bf;
y_MR = X * pa2;
NESE_MR = get_NMSE(y_MR, y0);    
fprintf('NMSE_MR = %.2f dB\n', NESE_MR);
%% SBL
[y_SBL,gain_rmsSBL] = get_SBL(K,M,X,H_x,x,y_bf,3);
NESE_SBL = get_NMSE(y_SBL, y0);
fprintf('NMSE_SBL = %.2f dB\n', NESE_SBL);