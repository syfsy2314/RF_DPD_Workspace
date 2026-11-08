function diag_real_coverage()
% 真实数据对照：LD-PLS 偶数阶基函数对
%   y0  = PA 输出（含 AM/AM + AM/PM）
%   y_f = PA 输出再过 OMUX（额外含频率相关相位）
% 的覆盖能力差异，直接量化 OMUX 相位对 LD-PLS 的额外伤害。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
x = x(:); y0 = y0(:);
K = 7; M = 5;

% OMUX 全通（与 verify_routeB 一致）
r = 0.8; theta = 0.6*pi;
b = [r^2, -2*r*cos(theta), 1]; a = [1, -2*r*cos(theta), r^2];
y_f = filter(b, a, y0);

% 对齐：分别把 y0、y_f 对齐到 x（互相关整数时延 + 复常数）
y0a = align_to(x, y0);
yfa = align_to(x, y_f);

% 偶数阶 BLMP 基函数（同相 + 记忆）
H = geth_even_1d(x, K, M);

y_hat_y0 = H * ((H'*H) \ (H'*y0a));
y_hat_yf = H * ((H'*H) \ (H'*yfa));

fprintf('--- 真实数据：LD-PLS 偶数阶覆盖能力 ---\n');
fprintf('拟合 y0  (PA, 含 AM/PM)      : NMSE = %6.2f dB\n', get_NMSE(y_hat_y0, y0a));
fprintf('拟合 y_f (PA + OMUX 相位)    : NMSE = %6.2f dB\n', get_NMSE(y_hat_yf, yfa));
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
