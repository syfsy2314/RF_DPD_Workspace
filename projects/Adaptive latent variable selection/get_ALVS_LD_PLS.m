function [y_ALVS, gain_rms, A] = get_ALVS_LD_PLS( ...
    K, M, X, x, h, y_bf, S, lamda, T_delta, kappa_G_max)
%GET_ALVS_LD_PLS 自适应潜变量选择的 LD-PLS 模型。
%   T_delta 控制新增潜变量的最小边际收益，kappa_G_max 控制潜变量
%   Gram 矩阵的最大条件数。触发停止条件的当前潜变量不会被保留。

H_x = geth_even_BLMP(x, K, M, h);

%% 分离线性分量
x_lin = H_x(:, 1);
X_nl = H_x(:, 2:end);
w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);
y_resid = y_bf - x_lin * w_lin;

proj_coeffs = (x_lin' * x_lin) \ (x_lin' * X_nl);
X0 = X_nl - x_lin * proj_coeffs;
Y0 = y_resid;

%% 标准化
x_mu = mean(X0, 1);
x_sigma = std(X0, 0, 1);
y_mu = mean(Y0);
y_sigma = std(Y0);
X_std = (X0 - x_mu) ./ x_sigma;
Y_std = (Y0 - y_mu) ./ y_sigma;

%% 自适应提取潜变量
Kmax = round(S * size(H_x, 2));
W = zeros(size(X_std, 2), Kmax);
P = zeros(size(X_std, 2), Kmax);
Z_score = zeros(size(X_std, 1), Kmax);
X_res = X_std;
Y_res = Y_std;
E_prev = real(Y_res' * Y_res);
A = 0;

for k = 1:Kmax
    w = X_res' * Y_res;
    w = w / norm(w);
    t = X_res * w;
    p = X_res' * t / (t' * t);
    q = Y_res' * t / (t' * t);

    X_res_new = X_res - t * p';
    Y_res_new = Y_res - t * q;
    E_new = real(Y_res_new' * Y_res_new);
    Delta_k = (E_prev - E_new) / E_prev;

    Z_try = [Z_score(:, 1:k-1), t];
    kappa_G_k = cond(Z_try' * Z_try);

        if kappa_G_k >= kappa_G_max && Delta_k <= T_delta
            break;
        end

    W(:, k) = w;
    P(:, k) = p;
    Z_score(:, k) = t;
    X_res = X_res_new;
    Y_res = Y_res_new;
    E_prev = E_new;
    A = k;
end
fprintf('选择的潜变量个数为 %d\n',A);
%% 潜变量空间岭回归
W = W(:, 1:A);
P = P(:, 1:A);
R = W / (P' * W);
Z = X_std * R;
C = (Z' * Z + lamda * eye(A)) \ (Z' * Y_std);
B_std = R * C;
B = B_std ./ x_sigma.' * y_sigma;

%% 使用全带 BLMP 基函数恢复 PA 输出
proj_coeff = (X(:, 1)' * X(:, 1)) \ ...
    (X(:, 1)' * X(:, 2:end));
X_nl_full = X(:, 2:end) - X(:, 1) * proj_coeff;
y_ALVS = X_nl_full * B + X(:, 1) * w_lin;
gain_rms = rms(y_ALVS) / rms(x);
end
