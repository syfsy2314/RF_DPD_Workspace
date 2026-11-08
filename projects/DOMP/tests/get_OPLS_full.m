function [Y_full, B_full] = get_OPLS_full(K, M, X, x, h, y_bf, lamda, O_max)
% get_OPLS_full  LD-PLS 增量扫描版: 潜变量一次提取, 返回全部 O 的重建与系数
%
%   [Y_full, B_full] = get_OPLS_full(K, M, X, x, h, y_bf, lamda, O_max)
%
%   与 get_OPLS 逐 O 调用的结果严格一致(PLS 潜变量前缀嵌套:
%   前 O 个潜变量与 O 的取值无关), 但潜变量提取只做一次,
%   复杂度从 O(T^2) 次完整调用降为 1 次提取 + T 次轻量岭求解。
%
%   输入:
%       K, M   - MP 阶数与记忆深度
%       X      - 全带基 (geth_even_1d), 用于重建
%       x      - 原始输入信号
%       h      - 反馈 FIR (geth_even_BLMP 用)
%       y_bf   - 带限输出 (辨识目标)
%       lamda  - 岭正则化参数(标量)
%       O_max  - 可选, 扫描的潜变量个数上限(默认 T-1 = 全部非线性列)
%   输出:
%       Y_full - N×O_max 复矩阵, 第 O 列 = 用前 O 个潜变量(+线性项)的全带重建
%       B_full - (T-1)×O_max, 第 O 列 = 前 O 个潜变量对应的非线性系数
%                (B_full 的行与 X(:,2:end) 列对应; 线性项系数 w_lin 未含)
%   其中 T = K*(M+1), O = round(S*T) 中的 T 即此值

if nargin < 8 || isempty(O_max)
    O_max = K*(M+1) - 1;
end

H_x = geth_even_BLMP(x, K, M, h);
Len = 1;
x_lin = H_x(:, Len);
X_nl_raw = H_x(:, (Len+1):end);

%% 线性剥离(与 get_OPLS 相同)
w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);
y_resid = y_bf - x_lin * w_lin;      % 剥离用带限线性输出(仅辨识域)

proj_coeffs = (x_lin' * x_lin) \ (x_lin' * X_nl_raw);
X_nl_clean = X_nl_raw - x_lin * proj_coeffs;

proj_coeff = (X(:,Len)' * X(:,Len)) \ (X(:,Len)' * X(:,(Len+1):end));
X_res = X(:, (Len+1):end) - X(:,Len) * proj_coeff;

O_max = size(X_nl_clean, 2);

%% PLS 潜变量提取(与 pls_B 相同流程, 一次提取到 O_max)
Xs = (X_nl_clean - mean(X_nl_clean))./std(X_nl_clean);
Ys = (y_resid - mean(y_resid))./std(y_resid);
W = zeros(size(Xs,2), O_max);
P = zeros(size(Xs,2), O_max);
Q = zeros(1, O_max);
X_resid = Xs;
Y_resid = Ys;
for k = 1:O_max
    w = X_resid' * Y_resid;
    w = w / norm(w);
    t = X_resid * w;
    p = X_resid' * t / (t'*t);
    q = Y_resid' * t / (t'*t);
    X_resid = X_resid - t * p';
    Y_resid = Y_resid - t * q;
    W(:,k) = w;
    P(:,k) = p;
    Q(k) = q;
end

%% 逐 O 增量求解(只用前 O 列, 岭回归)与全带重建
sx = std(X_nl_clean)';
sy = std(y_resid);
Y_full = zeros(length(y_bf), O_max);
B_full = zeros(O_max, O_max);
for O = 1:O_max
    R = W(:,1:O) * inv(P(:,1:O)' * W(:,1:O));
    Tm = Xs * R;
    C = inv(Tm'*Tm + lamda*eye(O)) * Tm' * Ys;
    B = R * C;
    B = B ./ sx * sy;                % 还原尺度(与 pls_B 一致)
    B_full(:, O) = B;
    % 重建线性项必须用【全带】列 X(:,1) 乘 w_lin(与 get_OPLS 第33行一致),
    % 带限列只用于辨识域的剥离, 不能用于重建
    Y_full(:, O) = X_res * B + X(:,Len) * w_lin;
end
end
