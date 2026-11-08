function diag_cyclo()
% 估计实测信号 y0 的循环频率（符号率），排除平稳成分后找峰值。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d);
N = numel(x);

% 频谱（看带宽和信号性质）
[p, f] = pwelch(y0, hann(1024), 512, 1024, 1);
[~, im] = max(p);
fprintf('y0 频谱峰值频率: %.4f (归一化), 能量带宽约 ±%.3f\n', f(im), f(im));

tau_max = min(512, N-1);
E = zeros(N,1);
for tau = 1:tau_max
    z = y0(tau+1:end) .* conj(y0(1:end-tau));
    E = E + abs(fft(z, N)).^2;
end

% 排除平稳成分：alpha≈0 (bin 1..16) 及镜像 alpha≈1 (bin N-15..N)
E([1:16, N-15:N]) = 0;

[sorted, idx] = sort(E, 'descend');
fprintf('排除平稳后 y0 循环频率峰值 alpha: %.5f %.5f %.5f %.5f %.5f\n', ...
    idx(1:5)/N);
fprintf('对应相对能量: %.3g %.3g %.3g %.3g %.3g\n', sorted(1:5)/max(sorted));
end
