function verify_routeB()
% 验证路B核心假设：PA 主幅度失真、OMUX 主相位失真，
% 即 频域相位差 angle(Y_f)-angle(X) 是否 ≈ OMUX 的真实相位。
% OMUX 用一个已知全通 IIR（r,theta）模拟相位失真。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
N = numel(x);

% OMUX = 全通 IIR（真值）
r = 0.8; theta = 0.6*pi;
b = [r^2, -2*r*cos(theta), 1];
a = [1, -2*r*cos(theta), r^2];
y_f = filter(b, a, y0);

% 简单对齐：互相关整数时延 + 常数相位补偿（不补偿频率相关群时延）
Hf = fft(x) .* conj(fft(y_f));
[~, ind] = max(abs(ifft(Hf)));
y_aligned = circshift(y_f, -(ind-1));
c = (x' * y_aligned) / (x' * x);
y_aligned = y_aligned * conj(c) / abs(c);

X = fftshift(fft(x));
Y = fftshift(fft(y_aligned));

% OMUX 真实相位
H1 = fftshift(freqz(b, a, N, 'whole'));
phi_true = angle(H1);

% 频域相位差
phi_est = angle(Y) - angle(X);

f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;      % 带内（信号带宽 50MHz/1GHz ≈ 0.05，取一半）
in = find(mask);

% 用"相位差"而非绝对相位来消除常数相位残留：dphi = phi(f)-phi(f-df)
df = 8;   % 频率差分步长（bin）
dphi_est  = phi_est(df+1:end) - phi_est(1:end-df);
dphi_true = phi_true(df+1:end) - phi_true(1:end-df);
ok = mask(df+1:end) & mask(1:end-df);
e_diff = mean(abs(wrapToPi(dphi_est(ok) - dphi_true(ok))));

e_abs = mean(abs(wrapToPi(phi_est(in) - phi_true(in))));

fprintf('绝对相位差误差(带内): %.4f rad\n', e_abs);
fprintf('相位差分误差(带内):   %.4f rad\n', e_diff);
end
