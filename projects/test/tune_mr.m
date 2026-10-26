%MR 超参数调优（H2 阻带 1e-6 配置）
%对图正则化权重 a1(系数图) 与 a2(样本图) 做二维对数网格扫描，
%判据为全带重建 NMSE(get_NMSE(y_MR, y0))，与 LD-PLS/DOMP 对比实验同口径
%图结构 Lb/Lg 只依赖 H_x，与 a1/a2 无关，循环外构建一次
%%
clear;
close all;
fs = 1000;
BW = 60;                 % 采集带宽（与 DOMP/test 对比实验一致，可改）
nfft = 1024;
window = hann(nfft);
%% 滤波器：H1 模拟域 80 dB（保持实际指标），H2 数字域 1e-6（深阻带）
Dpass1 = (10^(0.5/20) - 1) / (10^(0.5/20) + 1);
Dstop1 = 1e-4;
Dpass2 = 0.0057501127785;
Dstop2 = 1e-6;
dens = 20;
[N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(fs/2), [1 0], [Dpass1, Dstop1]);
H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
[N2, Fo, Ao, W] = firpmord([BW, BW+3]/(fs/2), [1 0], [Dpass2, Dstop2]);
H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
h = H2_FIR.Numerator;
%% 数据与观测链路
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
y_f = H1_FIR(y0);
[y_f, x, ~] = align_circular_signal(y_f, x);
y_bf = H2_FIR(y_f);
%% 基函数
K = 7; M = 5;
H_x = geth_even_BLMP(x, K, M, h);
X = geth_even_1d(x, K, M);
%% MR 图结构（只构建一次，不随 a1/a2 变化）
Kb = 2;  nb = 0.2;
Lb = get_Lb(H_x, Kb, nb);
Kg = 20; ng = 0.2; g = 1000;
G = H_x(1:g, :);
Lg = get_Lg(G, Kg, ng);
G_Lg_G = G' * Lg * G;                % 预计算样本图项
HxHx = H_x' * H_x;
Hxy  = H_x' * y_bf;
%% 基准（换 1e-6 前的 MR 原参数）
a1_ref = 1e-4; a2_ref = 1e-5;
pa2_ref = inv(HxHx + a1_ref*Lb + a2_ref*G_Lg_G) * Hxy;
nmse_ref = get_NMSE(X*pa2_ref, y0);
%% LD-PLS 基准（固定规则 S=0.3, lamda=10）
[y_LDPLS, ~] = get_OPLS(K, M, X, x, h, y_bf, 0.3, 10);
nmse_ldpls = get_NMSE(y_LDPLS, y0);
%% 二维对数网格扫描
a1_list = 10.^(-10:1:2);
a2_list = 10.^(-10:1:2);
n1 = numel(a1_list); n2 = numel(a2_list);
nmse_grid = zeros(n1, n2);
warning('off', 'MATLAB:nearlySingularMatrix');   % a1,a2 双双取到 1e-10 的角落
for i1 = 1:n1                                    % 等价于无正则化 LS，奇异警告属预期
    for i2 = 1:n2
        pa2 = inv(HxHx + a1_list(i1)*Lb + a2_list(i2)*G_Lg_G) * Hxy;
        nmse_grid(i1, i2) = get_NMSE(X*pa2, y0);
    end
end
warning('on', 'MATLAB:nearlySingularMatrix');
[nmse_best, idx] = min(nmse_grid(:));
[i1_best, i2_best] = ind2sub(size(nmse_grid), idx);
fprintf('--- MR 调参结果 (H2 阻带 1e-6, BW=%d MHz) ---\n', BW);
fprintf('MR 原参数   a1=%.0e, a2=%.0e : %.2f dB\n', a1_ref, a2_ref, nmse_ref);
fprintf('MR 最优     a1=%.1e, a2=%.1e : %.2f dB\n', ...
    a1_list(i1_best), a2_list(i2_best), nmse_best);
fprintf('LD-PLS 基准 S=0.3, lamda=10 : %.2f dB\n', nmse_ldpls);
%% 最优参数下的频谱对比备用变量
pa2_best = inv(HxHx + a1_list(i1_best)*Lb + a2_list(i2_best)*G_Lg_G) * Hxy;
y_MR_best = X * pa2_best;
%% 绘图：NMSE 热图（着色值即 dB）
figure;
imagesc(log10(a1_list), log10(a2_list), nmse_grid);
axis xy;
xlabel('log_{10}(a1) 系数图权重');
ylabel('log_{10}(a2) 样本图权重');
title('MR 全带重建 NMSE (dB) 参数热图');
colorbar;
hold on;
plot(log10(a1_list(i1_best)), log10(a2_list(i2_best)), 'rp', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'r');
plot(log10(a1_ref), log10(a2_ref), 'wo', 'MarkerSize', 10, 'LineWidth', 1.5);
legend({'最优', '原参数'}, 'Location', 'northeastoutside');
