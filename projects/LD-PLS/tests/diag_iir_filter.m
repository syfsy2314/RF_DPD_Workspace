function diag_iir_filter()
% 用 IIR 模拟滤波器（幅度+相位频率响应，Chebyshev-I 低通）代替全通 IIR，
% 验证：1) LD-PLS 偶数阶基函数对它的覆盖能力；
%       2) 路B 假设：带内频域响应比值 Y(f)/X(f) 是否 ≈ 滤波器 H_filt(f)。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
x = x(:); y0 = y0(:);
N = numel(x);
fs = 1000;                       % MHz，信号带宽 50MHz（±25MHz）

% --- IIR 模拟滤波器（低通，通带 35MHz / 阻带 55MHz，采样 1GHz）---
Wp = 35/(fs/2); Ws = 55/(fs/2); Rp = 0.5; Rs = 50;
[n, Wn] = cheb1ord(Wp, Ws, Rp, Rs);
[b, a] = cheby1(n, Rp, Wn);

% 滤波器频率响应（带内幅度+相位）
Hf = fftshift(freqz(b, a, N, 'whole'));
f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;           % 带内 ±30MHz
mag_in = abs(Hf(mask));
phi_in = unwrap(angle(Hf(mask)));
% 去线性（群时延）后的相位波动
P = polyfit(f(mask), phi_in, 1);
phi_dev = phi_in - polyval(P, f(mask));

fprintf('IIR 滤波器阶数 n=%d\n', n);
fprintf('带内幅度波动: %.3f dB (min %.3f, max %.3f)\n', ...
    20*log10(max(mag_in)/min(mag_in)), min(mag_in), max(mag_in));
fprintf('带内相位波动(去群时延): RMS %.2f deg, max %.2f deg\n', ...
    rms(phi_dev)*180/pi, max(abs(phi_dev))*180/pi);

% 对齐
y0a = align_to(x, y0);
y_f = filter(b, a, y0a);
y_fa = align_to(x, y_f);

% --- 1) LD-PLS 偶数阶基函数覆盖能力 ---
K = 7; M = 5;
H = geth_even_1d(x, K, M);
y_hat0 = H * ((H'*H) \ (H'*y0a));
y_hatf = H * ((H'*H) \ (H'*y_fa));
fprintf('\nLD-PLS 拟合 y0  (PA):      NMSE = %6.2f dB\n', get_NMSE(y_hat0, y0a));
fprintf('LD-PLS 拟合 y_f (PA+IIR):  NMSE = %6.2f dB\n', get_NMSE(y_hatf, y_fa));

% --- 2) 路B：带内频域响应比值估计 ---
X = fftshift(fft(x));
Y = fftshift(fft(y_fa));
Hest = Y ./ X;                    % 频域比值
mag_est = abs(Hest(mask));
phi_est = unwrap(angle(Hest(mask)));

% 幅度：各自去均值后比较波动形状（消除整体增益差异）
mag_true_n = mag_in / mean(mag_in);
mag_est_n  = mag_est / mean(mag_est);
e_mag = mean(abs(mag_est_n - mag_true_n));

% 相位：去线性相位后比较
Pe = polyfit(f(mask), phi_est, 1);
phi_est_dev = phi_est - polyval(Pe, f(mask));
e_phi = mean(abs(wrapToPi(phi_est_dev - phi_dev)));

fprintf('\n路B 频域比值估计误差（带内）:\n');
fprintf('  归一化幅度误差: %.4f\n', e_mag);
fprintf('  相位误差(去群时延): %.3f deg\n', e_phi*180/pi);
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
