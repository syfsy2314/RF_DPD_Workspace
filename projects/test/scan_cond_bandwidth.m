%用于测试随着滤波器带宽的降低，条件数与线性剥离效果的变化（带宽扫描版）
%输出：cond(H_x)、cond(H_n) 两条条件数曲线，以及线性剥离前后第一潜变量与线性列的相关性
%%
clear;
close all;
fs = 1000;  % 采样率（MHz），与 BW 单位一致
%% 导入数据（与带宽无关，只加载一次）
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
%% 扫描设置
BW_list = 100:-5:30;            % 采集带宽扫描（MHz），从大到小
nBW = numel(BW_list);
kappa_x = zeros(1, nBW);         % 原始基条件数（LS 问题口径）
kappa_n = zeros(1, nBW);         % 列归一化基条件数（纯共线性口径）
r_raw   = zeros(1, nBW);         % 不剥离线性：第一潜变量与线性列的相关系数
r_det   = zeros(1, nBW);         % 剥离线性后：第一潜变量与线性列的相关系数
sig_max = zeros(1, nBW);         % 最大奇异值（主导列能量）
sig_min = zeros(1, nBW);         % 最小奇异值（反映共线性）
K = 7;
M = 5;
Dstop1 = 0.0001;
Dstop2 = 0.000001;
dens  = 20;
Dpass1 = (10^(0.5/20) - 1) / (10^(0.5/20) + 1);      % H1 通带波纹 0.5 dB
Dpass2 = 0.0057501127785;                            % H2 通带波纹（原脚本取值）
%% 带宽扫描
for iBW = 1:nBW
    BW = BW_list(iBW);
    % 第一个 FIR 滤波器（OMUX 等效，通带随 BW 缩放）
    [N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(fs/2), [1 0], [Dpass1, Dstop1]);
    H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
    % 第二个 FIR 滤波器（反馈覆盖）
    [N2, Fo, Ao, W] = firpmord([BW, BW+3]/(fs/2), [1 0], [Dpass2, Dstop2]);
    H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
    h2 = H2_FIR.Numerator;
    % 观测链路：把观测 y_f 对齐到输入 x（x 作参考不变，x_al 即 x 原样）
    y_f = H1_FIR(y0);
    [y_f, x_al, ~] = align_circular_signal(y_f, x);
    % 带限观测目标（与 LD-PLS 建模口径一致）
    y_bf = H2_FIR(y_f);
    % 带限基函数及其条件数
    H_x = geth_even_BLMP(x_al, K, M, h2);
    Hn = H_x ./ vecnorm(H_x);                        % 列归一化
    s1 = svd(H_x, 'econ');
    s2 = svd(Hn,  'econ');
    kappa_x(iBW) = s1(1) / max(s1(end), eps);   % eps 下限防止奇异值为 0 时出现 Inf
    kappa_n(iBW) = s2(1) / max(s2(end), eps);
    sig_max(iBW) = s1(1);                            % 最大奇异值
    sig_min(iBW) = s1(end);                          % 最小奇异值
    % 量化线性剥离效果：不剥离 vs 剥离线性下，第一潜变量与线性列的相关系数
    x_lin = H_x(:, 1);                                 % 线性基（BLMP 第 1 列）
    X_nl_raw = H_x(:, 2:end);                          % 非线性基（其余列）
    w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);        % 线性最小二乘增益
    y_resid = y_bf - x_lin * w_lin;                    % 剥离线性后的目标
    proj = (x_lin' * x_lin) \ (x_lin' * X_nl_raw);     % 各非线性列在 x_lin 上的投影系数
    X_nl_clean = X_nl_raw - x_lin * proj;              % 剥离线性后的非线性基
    r_tmp = pls_dir_corr(H_x, y_bf, x_lin);            % 不剥离：完整基做 PLS
    r_raw(iBW) = r_tmp(1);                             % 第一潜变量与线性列相关
    r_tmp = pls_dir_corr(X_nl_clean, y_resid, x_lin);  % 剥离后：残差基做 PLS
    r_det(iBW) = r_tmp(1);
    fprintf('BW = %3d MHz: cond(H_x) = %.3g, cond(H_n) = %.3g, r_raw = %.3f, r_det = %.3f\n', ...
        BW, kappa_x(iBW), kappa_n(iBW), r_raw(iBW), r_det(iBW));
end
%% 绘图（纵轴为 log10(条件数)，线性坐标）
figure;
plot(BW_list, log10(kappa_x), 'o-', 'LineWidth', 1.2); hold on;
plot(BW_list, log10(kappa_n), 's-', 'LineWidth', 1.2);
xlabel('采集带宽 (MHz)');
ylabel('log_{10}(条件数)');
set(gca, 'XDir', 'reverse');         % 横轴从大到小显示（100 -> 30）
legend('\kappa(H_x) 原始基', '\kappa(H_n) 列归一化基', 'Location', 'best');
title('BLMP 基函数矩阵条件数随采集带宽的变化');
grid on;
%% 线性剥离效果：不剥离 vs 剥离下第一潜变量与线性列的相关性
figure;
plot(BW_list, r_raw, 'o-', 'LineWidth', 1.2); hold on;
plot(BW_list, r_det, 's-', 'LineWidth', 1.2);
xlabel('采集带宽 (MHz)');
ylabel('|cos(第一潜变量, x_{lin})|');
set(gca, 'XDir', 'reverse');         % 横轴从大到小显示（100 -> 30）
legend('不剥离（第一潜变量≈线性）', '剥离后（第一潜变量对准非线性）', 'Location', 'best');
title('线性剥离对第一潜变量方向的影响随带宽的变化');
grid on;
%% 奇异值分解：定位条件数下降由谁主导
figure;
semilogy(BW_list, sig_max, 'o-', 'LineWidth', 1.2); hold on;
semilogy(BW_list, sig_min, 's-', 'LineWidth', 1.2);
xlabel('采集带宽 (MHz)');
ylabel('奇异值（对数坐标）');
set(gca, 'XDir', 'reverse');         % 横轴从大到小显示（100 -> 30）
legend('\sigma_{max}（主导列能量）', '\sigma_{min}（共线性）', 'Location', 'best');
title('最大/最小奇异值随带宽的变化：判断条件数下降的机理');
grid on;

function r = pls_dir_corr(X0, Y0, x_lin)
% 跑完整 NIPALS，返回每个潜变量得分向量 t_k 与线性列 x_lin 的余弦相似度
X = (X0 - mean(X0)) ./ std(X0);
Y = (Y0 - mean(Y0)) ./ std(Y0);
xl = (x_lin - mean(x_lin)) / std(x_lin);
X_res = X; Y_res = Y;
K = size(X, 2);
r = zeros(1, K);
for k = 1:K
    w = X_res' * Y_res;
    w = w / norm(w);
    t = X_res * w;
    r(k) = abs(t' * xl) / (norm(t) * norm(xl));
    p = X_res' * t / (t' * t);
    q = Y_res' * t / (t' * t);
    X_res = X_res - t * p';
    Y_res = Y_res - t * q;
end
end
