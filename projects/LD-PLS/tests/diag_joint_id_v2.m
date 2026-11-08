function diag_joint_id_v2()
% 路B 理论实现：过参数化 Hammerstein 辨识 + 交叉关系法分离。
% 系统: x -> PA(记忆多项式) -> y0 -> H(椭圆 IIR, 未知) -> y_f
% 仅用 (x, y_f)：先过参数化线性回归，再用多通道盲解卷积分离 H 与 PA。

% ---- 数据 ----
d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d); y0 = y0(:);
N = numel(x);

% ---- H 真值：椭圆 IIR（辨识时视为未知）----
fs = 1000;
Wp = 35/(fs/2); Ws = 55/(fs/2); Rp = 0.5; Rs = 50;
[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b0, a0] = ellip(n, Rp, Rs, Wn);
fprintf('椭圆 IIR 阶数 n=%d\n', n);

% ---- y_f（带限）----
y_f = align_to(x, filter(b0, a0, y0));

% ---- 参数 ----
K = 7;          % 非线性阶数
M = 5;          % PA 记忆深度
P = 60;         % H 的 FIR 近似长度
L = P - 1 + M;  % 总延迟
Ma = M + 1;     % a_k 长度
Ld = L + 1;     % d_k 长度 = P + M

% ---- 1) 过参数化线性回归（只依赖 x，稳健）----
Psi = zeros(N, Ld*K);
col = 0;
for l = 0:L
    xl = [zeros(l,1); x(1:end-l)];
    axl = abs(xl);
    for k = 1:K
        col = col + 1;
        Psi(:, col) = xl .* axl.^(k-1);
    end
end
dvec = (Psi'*Psi) \ (Psi'*y_f);
D = reshape(dvec, Ld, K);           % D(:,k) = d_k = h ⊛ a_k

% ---- 2) 交叉关系法分离 ----
% 交叉关系: a_1 ⊛ d_k = a_k ⊛ d_1  =>  T(d_1)a_k - T(d_k)a_1 = 0
nb_rows = (Ma + Ld - 1) * (K-1);
A_mat = zeros(nb_rows, K*Ma);
for k = 2:K
    r0 = (k-2)*(Ma + Ld - 1);
    T1 = convmat(D(:,1), Ma);          % T(d_1)
    Tk = convmat(D(:,k), Ma);          % T(d_k)
    A_mat(r0+1:r0+Ma+Ld-1, 1:Ma)        = -Tk;
    A_mat(r0+1:r0+Ma+Ld-1, (k-1)*Ma+1:k*Ma) = T1;
end
[~, ~, V] = svd(A_mat, 'econ');
avec = V(:, end);
avec = avec / avec(1);                % 归一化 a_{1,0} = 1（消除尺度模糊）
A_res = reshape(avec, Ma, K);         % A_res(:,k) = a_k

% ---- 3) 恢复 h（FIR 近似 H）----
a1 = A_res(:, 1);
h = convmat(a1, P) \ D(:, 1);         % d_1 = h ⊛ a_1

% ---- 验证 1：PA 重构 ----
% geth_even_1d 列序: m 外层 k 内层，与 reshape(A_res', [], 1) 一致
Phi = geth_even_1d(x, K, M);
a_pa = reshape(A_res', [], 1);
y0_hat = Phi * a_pa;
fprintf('PA 拟合 y0      : NMSE = %6.2f dB\n', get_NMSE(y0_hat, y0));

% ---- 验证 2：H 频率响应（归一化带内比较）----
H_true = fftshift(freqz(b0, a0, N, 'whole'));
H_hat  = fftshift(freqz(h,  1,  N, 'whole'));
f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;
H_true = H_true / H_true(find(mask,1));
H_hat  = H_hat  / H_hat(find(mask,1));
e_mag = mean(abs(abs(H_hat(mask)) - abs(H_true(mask))));
e_phi = mean(abs(wrapToPi(angle(H_hat(mask)) - angle(H_true(mask)))));
fprintf('H 带内幅度误差 : %.4f\n', e_mag);
fprintf('H 带内相位误差 : %.3f deg\n', e_phi*180/pi);
end

function T = convmat(d, len)
% 返回 T 使 T*x = conv(d, x)，x 长度为 len。
Ld = length(d);
T = zeros(Ld + len - 1, len);
for i = 1:len
    T(i:i+Ld-1, i) = d;
end
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
