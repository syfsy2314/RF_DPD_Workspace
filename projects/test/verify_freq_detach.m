function verify_freq_detach()
% 频域比值去卷积 + PA 拟合：验证"半盲分离"的核心可行性。
% 1) H_est = Y_f/X（带内，非参数）
% 2) 去卷积恢复 y0_est = ifft(Y_f/H_est)（带内，带外正则化）
% 3) 拟合 PA（记忆多项式）到 y0_est，比较 NMSE(y0_hat, y0)

clear;
d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d); y0 = y0(:);
N = numel(x);

% 椭圆 IIR 滤波器（低通，覆盖双频段信号）
[b0,a0]=ellip(6,0.5,50,0.25);
y_f = align_to(x, filter(b0,a0,y0));

K=7; M=5;
Phi = geth_even_1d(x, K, M);

f=(-N/2:N/2-1).'/N;
X=fftshift(fft(x));
Y=fftshift(fft(y_f));

% 频域比值估计 H（带内）
B = 0.25;
mask = abs(f) <= B;
H_est = Y ./ X;
H_est(~mask) = 0;                    % 带外置 0

% 去卷积恢复 y0（带内）
Y0_freq = zeros(N,1);
Y0_freq(mask) = Y(mask) ./ H_est(mask);
y0_est = ifft(ifftshift(Y0_freq));

% 拟合 PA
a = (Phi'*Phi) \ (Phi'*y0_est);
y0_hat = Phi * a;

fprintf('频域去卷积恢复 y0: NMSE = %.2f dB\n', get_NMSE(y0_est, y0));
fprintf('PA 拟合 y0_hat vs y0: NMSE = %.2f dB\n', get_NMSE(y0_hat, y0));
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
