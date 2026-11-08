function diag_sym()
% 精确检测信号的符号率 / 循环前缀结构（物理先验）。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
N = numel(x);

% --- 1) 瞬时功率 |x|^2 的循环频率（符号率/脉冲成形相关）---
p = abs(x).^2;
P = abs(fft(p, N));
P = P(1:floor(N/2));
P(1:40) = 0;                          % 排除直流附近
[~, im] = max(P);
fprintf('|x|^2 循环频率峰值 alpha=%.5f  (相对能量 %.3f)\n', im/N, P(im)/sum(P));

% --- 2) 循环自相关函数 |R_x^alpha(0)|，扫 alpha ---
alphas = (0.001:0.0005:0.2);          % 扫描符号率范围
R = zeros(size(alphas));
n = (0:N-1).';
for i = 1:numel(alphas)
    a = alphas(i);
    R(i) = abs(mean(x .* conj(x) .* exp(-1j*2*pi*a*n)));
end
R = R / max(R);
% 打印前几个局部峰
[pk, loc] = findpeaks(R, 'MinPeakHeight', 0.3);
fprintf('循环自相关峰值 alpha: ');
if isempty(pk)
    fprintf('无显著峰\n');
else
    for i = 1:numel(loc)
        fprintf('%.4f(%.2f) ', alphas(loc(i)), pk(i));
    end
    fprintf('\n');
end

% --- 3) 大范围自相关找 CP 周期（全 N，FFT 加速）---
R2 = ifft(abs(fft(x, 2*N)).^2);
R2 = abs(R2(1:N)) / abs(R2(1));
% 排除 tau < 10 的近邻，找最大次峰
R2(1:10) = 0;
[~, tau] = max(R2);
fprintf('自相关最大次峰 tau=%d (幅度 %.4f)\n', tau, R2(tau));
end
