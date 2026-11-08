function diag_pa_phase()
% 分解 PA 频域相位差：线性部分（分数时延/常数）与非线性部分（AM/PM）。
% 目的：判断 verify_routeB_control 中 PA 的 0.64 rad 尖峰是时延残留还是真非线性。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
N = numel(x);

% 对齐（与 verify_routeB 一致）
Hf = fft(x) .* conj(fft(y0));
[~, ind] = max(abs(ifft(Hf)));
y0a = circshift(y0, -(ind-1));
c = (x' * y0a) / (x' * x);
y0a = y0a * conj(c) / abs(c);

X = fftshift(fft(x));
Y0 = fftshift(fft(y0a));
phi_pa = angle(Y0) - angle(X);

f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;
in = find(mask);
fi = f(in);
pi_in = unwrap(phi_pa(in));

% 线性拟合 phi = a*f + b，残差即非线性 AM/PM 相位
P = polyfit(fi, pi_in, 1);
lin = polyval(P, fi);
resid = pi_in - lin;

fprintf('--- PA 相位分解（带内）---\n');
fprintf('线性斜率: %.4f rad/归一化频率 (对应分数时延约 %.2f 采样)\n', ...
    P(1), P(1)/(2*pi));
fprintf('线性相位带内幅度范围: %.4f rad\n', max(lin)-min(lin));
fprintf('非线性残差 RMS=%.4f rad, max|.|=%.4f rad\n', ...
    rms(resid), max(abs(resid)));
fprintf('非线性 / 线性 (RMS 比): %.3f\n', rms(resid) / rms(lin - mean(lin)));
end
