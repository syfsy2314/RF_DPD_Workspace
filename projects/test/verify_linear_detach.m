%% 验证线性剥离的必要性
% 证据1：潜变量方向被线性霸占（不剥离时第一潜变量与线性列相关系数≈1）
% 证据2：同参数下，剥离 vs 不剥离的 NMSE 收益
clear; close all;
fs = 1000;
BW = 50;
K = 7; M = 5;
S = 0.3; lamda = 10;
T = K * (M + 1);
O = round(S * T);
%% 加载数据
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
%% 滤波器（与 main.m 一致）
Fs = 1000;
Fpass = BW + 5; Fstop = BW + 10;
Dpass = (10^(0.5/20) - 1) / (10^(0.5/20) + 1); Dstop = 0.0001; dens = 20;
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
H1_FIR = dsp.FIRFilter('Numerator', firpm(N, Fo, Ao, W, {dens}));
Fpass = BW; Fstop = BW + 3;
Dpass = 0.0057501127785; Dstop = 1e-6;
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
H2_FIR = dsp.FIRFilter('Numerator', firpm(N, Fo, Ao, W, {dens}));
h2 = H2_FIR.Numerator;
%% 观测链路
y_f = H1_FIR(y0);
[y_f, x_al, ~] = align_circular_signal(y_f, x);
y_bf = H2_FIR(y_f);
H_x = geth_even_BLMP(x_al, K, M, h2);
X   = geth_even_1d(x_al, K, M);
%% 线性剥离
x_lin = H_x(:, 1);
X_nl_raw = H_x(:, 2:end);
w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);
y_resid = y_bf - x_lin * w_lin;
proj = (x_lin' * x_lin) \ (x_lin' * X_nl_raw);
X_nl_clean = X_nl_raw - x_lin * proj;
%% 证据1：潜变量方向与线性列的相关性
r_raw = pls_dir_corr(H_x, y_bf, x_lin);          % 不剥离
r_det = pls_dir_corr(X_nl_clean, y_resid, x_lin); % 剥离后
fprintf('不剥离:  第1潜变量与线性列相关 = %.4f\n', r_raw(1));
fprintf('剥离后:  第1潜变量与线性列相关 = %.4f\n', r_det(1));
figure;
plot(1:numel(r_raw), r_raw, 'o-', 'LineWidth', 1.2); hold on;
plot(1:numel(r_det), r_det, 's-', 'LineWidth', 1.2);
xlabel('潜变量序号 k');
ylabel('|cos(t_k, x_{lin})|');
legend('不剥离（第1个潜变量≈线性项）', '剥离后（潜变量对准非线性）', 'Location', 'best');
title('线性剥离把 PLS 潜变量从线性方向释放出来');
grid on;
%% 证据2：同参数 NMSE 对比（唯一差异是是否剥离）
B_nodt = pls_B(H_x, y_bf, O, lamda);        % 不剥离，同样岭回归
y_nodt = X * B_nodt;
NESE_nodt = get_NMSE(y_nodt, y0);
[y_dt, ~] = get_OPLS(K, M, X, x_al, h2, y_bf, S, lamda);
NESE_dt = get_NMSE(y_dt, y0);
fprintf('================================================\n');
fprintf('不剥离 PLS  : NMSE = %.2f dB\n', NESE_nodt);
fprintf('线性剥离 PLS : NMSE = %.2f dB\n', NESE_dt);
fprintf('================================================\n');

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
