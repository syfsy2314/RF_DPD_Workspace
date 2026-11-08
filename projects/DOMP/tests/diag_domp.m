%% 诊断: 为什么 DOMP 带内拟合好但全带重建差
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
%% 基函数
K = 7; M = 5;
H_x = geth_even_BLMP(x, K, M, h);
X   = geth_even_1d(x, K, M);
nrm = @(a) sprintf('%8.2f dB', get_NMSE(a, y0));
%% 1) 条件数
fprintf('cond(X''X)   = %.3g\n', cond(X'*X));
fprintf('cond(Hx''Hx) = %.3g\n', cond(H_x'*H_x));
%% 2) 全 LS 基准(无正则化)
a_ls = pinv(H_x)*y_bf;
fprintf('全LS(42项) 重建: %s\n', nrm(X*a_ls));
%% 3) 岭回归 lamda 扫描
fprintf('--- 岭回归扫描 ---\n');
sc = norm(H_x'*H_x)/size(H_x,2);   % Gram 尺度
for lam = [1e-6 1e-5 1e-4 1e-3 1e-2 0.1]*sc
    a_r = (H_x'*H_x + lam*eye(size(H_x,2))) \ (H_x'*y_bf);
    fprintf('lam=%.1e: 重建 %s\n', lam, nrm(X*a_r));
end
%% 4) DOMP 阶数扫描(固定口径, 不用 BIC)
Rmat_label = cell(1, K*(M+1));
[h_domp, s, nopt, h_full] = omp_domp(H_x, y_bf, Rmat_label, 'DOMP');
fprintf('BIC 选的 nopt = %d\n', nopt);
fprintf('--- DOMP 阶数扫描 ---\n');
for n = [5 8 12 17 22 30 42]
    fprintf('n=%2d: 重建 %s\n', n, nrm(X*h_full(:,n)));
end
%% 5) OMP 对照(不做字典正交化)
[h_omp, s_omp, nopt_omp, hfull_omp] = omp_domp(H_x, y_bf, Rmat_label, 'OMP');
fprintf('OMP: nopt=%d, 重建 %s\n', nopt_omp, nrm(X*h_omp));
%% 6) 选中的项(对照 geth_even_1d 列序: 外层m, 内层k)
fprintf('--- DOMP 选中的核(m, k) ---\n');
fprintf('%s\n', num2str([mod(s-1,K)+1, floor((s-1)/K)]));
