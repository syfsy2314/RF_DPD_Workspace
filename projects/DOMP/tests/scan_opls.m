%% LD-PLS (OPLS) 定阶扫描 —— 与 demo.m 中 DOMP 扫描同口径
% 扫描潜变量个数 O 与正则化 lamda, 判据均为全带重建 NMSE(对 y0)
% 注: 扫描判据用了 y0, 属"性能上界"口径; 与 DOMP 扫描对比时双方口径一致
clear; close all;
%% 滤波器(与 demo.m 相同)
BW = 60; Fs = 1000;
Rp_dB = 0.5; Dpass = (10^(Rp_dB/20)-1)/(10^(Rp_dB/20)+1); Dstop = 0.0001; dens = 20;
[N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(Fs/2), [1 0], [Dpass, Dstop]);
H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
[N2, Fo, Ao, W] = firpmord([BW, BW+3]/(Fs/2), [1 0], [Dpass, Dstop]);
H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
h = H2_FIR.Numerator;
%% 数据
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
y_f = H1_FIR(y0);
[y_f, x, ~] = align_circular_signal(y_f, x);
y_bf = H2_FIR(y_f);
%% 基函数与常量
K = 7; M = 5;
X = geth_even_1d(x, K, M);
T_all = K*(M+1);                         % 42, 其中非线性列 41 个
%% 1) O 扫描(lamda 固定)
lamda = 10;
nmse_O = zeros(1, T_all-1);              % O = 1..41
for O = 1:T_all-1
    y_o = get_OPLS(K, M, X, x, h, y_bf, O/T_all, lamda);
    nmse_O(O) = get_NMSE(y_o, y0);
end
[~, O_best] = min(nmse_O);
fprintf('--- O 扫描 (lamda=%g) ---\n', lamda);
fprintf('固定规则 O=%d: %.2f dB\n', round(0.3*T_all), nmse_O(round(0.3*T_all)));
fprintf('扫描最优 O=%d: %.2f dB\n', O_best, nmse_O(O_best));
%% 2) O × lamda 二维扫描
lam_list = [0.1 1 10 100 1000];
nmse_2d = zeros(T_all-1, numel(lam_list));
for il = 1:numel(lam_list)
    for O = 1:T_all-1
        y_o = get_OPLS(K, M, X, x, h, y_bf, O/T_all, lam_list(il));
        nmse_2d(O, il) = get_NMSE(y_o, y0);
    end
    fprintf('lamda=%g: 最优 O=%d, %.2f dB\n', lam_list(il), ...
        find(nmse_2d(:,il)==min(nmse_2d(:,il)),1), min(nmse_2d(:,il)));
end
[nmse_min, idx] = min(nmse_2d(:));
[O_gbest, il_gbest] = ind2sub(size(nmse_2d), idx);
fprintf('--- 二维全局最优: O=%d, lamda=%g -> %.2f dB ---\n', ...
    O_gbest, lam_list(il_gbest), nmse_min);
%% 3) 绘图
figure;
plot(1:T_all-1, nmse_2d, 'LineWidth', 1); hold on;
plot(O_gbest, nmse_min, 'ro', 'MarkerFaceColor', 'r');
xlabel('潜变量个数 O'); ylabel('重建 NMSE (dB)');
title('LD-PLS 定阶扫描(不同 lamda)');
legend(arrayfun(@(l) sprintf('lamda=%g', l), lam_list, 'UniformOutput', false), ...
    'Location', 'best');
grid on;
