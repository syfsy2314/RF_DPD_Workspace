function diag_omux_sensitivity()
% 敏感性扫描：OMUX 全通极点强度（r 逼近单位圆、theta 落在带内）
% 对 LD-PLS 偶数阶基函数覆盖能力的影响。
% 目的：找到"OMUX 相位失真多大时，LD-PLS 覆盖才真正失效、路B才值得做"的阈值。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
x = x(:); y0 = y0(:);
N = numel(x);
K = 7; M = 5;

y0a = align_to(x, y0);
H = geth_even_1d(x, K, M);

f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;

% 基线：拟合纯 PA
y_hat0 = H * ((H'*H) \ (H'*y0a));
nmse0 = get_NMSE(y_hat0, y0a);

fprintf('基线 拟合 y0 (PA): NMSE = %.2f dB\n\n', nmse0);
fprintf('  r     theta    OMUX相位RMS   OMUX相位max    拟合y_f NMSE    劣化\n');

% 把极点放在带内（theta 对应带内频率），扫 r 逼近单位圆
theta = 0.05*pi;              % f=theta/2pi = 0.025，带内中心附近
for r = [0.7, 0.85, 0.92, 0.96, 0.98, 0.995]
    b = [r^2, -2*r*cos(theta), 1];
    a = [1, -2*r*cos(theta), r^2];
    y_f = filter(b, a, y0a);
    y_fa = align_to(x, y_f);

    % OMUX 带内相位统计
    H1 = fftshift(freqz(b, a, N, 'whole'));
    phi = angle(H1); phi_in = phi(mask);
    phi_rms = rms(wrapToPi(phi_in));
    phi_max = max(abs(wrapToPi(phi_in)));

    y_hat = H * ((H'*H) \ (H'*y_fa));
    nmse = get_NMSE(y_hat, y_fa);

    fprintf('%5.3f  %5.3f   %8.2f deg   %8.2f deg    %7.2f dB    %+6.2f dB\n', ...
        r, theta, phi_rms*180/pi, phi_max*180/pi, nmse, nmse-nmse0);
end
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
