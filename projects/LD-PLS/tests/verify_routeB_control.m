function verify_routeB_control()
% 对照组：量化 PA 自身（y0，未过 OMUX）对频域相位的贡献。
% 若 PA 自身相位差 << OMUX 真实相位差，则坐实"PA 主幅度、OMUX 主相位"分离假设。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
N = numel(x);

% OMUX 真值相位（与 verify_routeB 一致的全通 IIR）
r = 0.8; theta = 0.6*pi;
b = [r^2, -2*r*cos(theta), 1];
a = [1, -2*r*cos(theta), r^2];
H1 = fftshift(freqz(b, a, N, 'whole'));
phi_omux = angle(H1);

% PA 自身：对齐 y0 到 x（互相关整数时延 + 常数相位补偿）
Hf = fft(x) .* conj(fft(y0));
[~, ind] = max(abs(ifft(Hf)));
y0a = circshift(y0, -(ind-1));
c = (x' * y0a) / (x' * x);
y0a = y0a * conj(c) / abs(c);

X = fftshift(fft(x));
Y0 = fftshift(fft(y0a));
phi_pa = angle(Y0) - angle(X);

f = (-N/2:N/2-1).'/N;
mask = abs(f) <= 0.03;      % 带内
in = find(mask);

phi_pa_in   = wrapToPi(phi_pa(in));
phi_omux_in = wrapToPi(phi_omux(in));

fprintf('--- 对照组：带内相位统计 ---\n');
fprintf('PA 自身相位差   RMS=%.4f rad, max|.|=%.4f rad\n', ...
    rms(phi_pa_in), max(abs(phi_pa_in)));
fprintf('OMUX 真值相位差 RMS=%.4f rad, max|.|=%.4f rad\n', ...
    rms(phi_omux_in), max(abs(phi_omux_in)));
fprintf('PA相位 / OMUX相位 (RMS 比): %.3f\n', ...
    rms(phi_pa_in) / rms(phi_omux_in));
end
