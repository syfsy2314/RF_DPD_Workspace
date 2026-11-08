function diag_joint_id()
% 联合辨识 PA(记忆多项式) + IIR 滤波器 H，仅用 (x, y_f)。
% 系统: x -> PA(记忆多项式) -> y0 -> IIR(B/A) -> y_f(带限)
% 方法: 交替迭代 —— 固定 H 逆滤波后拟合 PA；固定 PA 后用方程误差法辨识 IIR。

% ---- 数据 ----
d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d); y0 = y0(:);
N = numel(x);

% ---- IIR 真值（仿真用，辨识时视为未知）----
fs = 1000;
Wp = 35/(fs/2); Ws = 55/(fs/2); Rp = 0.5; Rs = 50;
[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b_true, a_true] = ellip(n, Rp, Rs, Wn);   % 椭圆滤波器
Nb = length(b_true)-1; Na = length(a_true)-1;
fprintf('椭圆 IIR 阶数 n=%d\n', n);

% ---- y_f（带限）----
y_f = align_to(x, filter(b_true, a_true, y0));

% ---- PA 偶数阶基函数 ----
K = 7; M = 5;
Phi = geth_even_1d(x, K, M);

% ---- 交替迭代辨识 ----
% 初始化 H：整体线性 IIR 辨识 (x -> y_f)
[b_hat, a_hat] = id_iir(x, y_f, Nb, Na);

nmse_hist = zeros(1, 30);
for iter = 1:30
    % 1) 固定 H：逆滤波恢复 y0，再拟合 PA
    y0_est = align_to(x, filter(a_hat, b_hat, y_f));   % H^{-1} * y_f
    a_pa   = (Phi'*Phi) \ (Phi'*y0_est);
    y0_hat = Phi * a_pa;

    % 2) 固定 PA：从 (y0_hat, y_f) 辨识 IIR
    [b_hat, a_hat] = id_iir(y0_hat, y_f, Nb, Na);

    % 收敛指标
    y_f_hat = filter(b_hat, a_hat, y0_hat);
    nmse_hist(iter) = get_NMSE(y_f_hat, y_f);
end

% ---- 评估 ----
fprintf('PA 拟合 y0   : NMSE = %6.2f dB\n', get_NMSE(y0_hat, y0));
fprintf('端到端 y_f   : NMSE = %6.2f dB\n', nmse_hist(end));

% H 频率响应比较（归一化 DC 增益，只看带内）
H_true = fftshift(freqz(b_true, a_true, N, 'whole'));
H_hat  = fftshift(freqz(b_hat,  a_hat,  N, 'whole'));
f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;
H_true = H_true / H_true(mask);   % 归一化（用带内中心 bin 的复增益）
H_hat  = H_hat  / H_hat(mask);
e_mag = mean(abs(abs(H_hat(mask)) - abs(H_true(mask))));
e_phi = mean(abs(wrapToPi(angle(H_hat(mask)) - angle(H_true(mask)))));
fprintf('H 带内幅度误差: %.4f\n', e_mag);
fprintf('H 带内相位误差: %.3f deg\n', e_phi*180/pi);
fprintf('迭代 NMSE 历史: '); fprintf('%.1f ', nmse_hist(1:5:end)); fprintf('dB\n');
end

function [b, a] = id_iir(u, y, nb, na)
% 方程误差法 IIR 辨识: y = B/A * u，参数化 a0 = 1。
N = length(u);
Z = zeros(N, na + nb + 1);
for i = 1:na
    Z(:, i) = -[zeros(i,1); y(1:end-i)];
end
for i = 0:nb
    Z(:, na+1+i) = [zeros(i,1); u(1:end-i)];
end
theta = Z \ y;
b = theta(na+1:end);
a = [1; theta(1:na)];
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
