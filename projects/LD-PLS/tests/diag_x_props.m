function diag_x_props()
% 确认 x 的性质：采样率、带宽、符号周期（循环频率）。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
N = numel(x);

% 频谱（归一化频率轴，默认采样率按 1 处理）
[p, f] = pwelch(x, hann(1024), 512, 1024, 'centered');
[~, im] = max(p);
fprintf('N=%d, 频谱峰值 f=%.4f (归一化)\n', N, f(im));

% 自相关（找循环前缀相关的符号周期峰值）
R = xcorr(x, x, min(round(N/2), 40000));
R = abs(R);
Rpos = R(round(length(R)/2)+2 : end);   % 正延迟 tau>=1
[sorted, idxs] = sort(Rpos, 'descend');
fprintf('自相关前 5 个峰值 tau: %d %d %d %d %d\n', idxs(1:5));
fprintf('峰值相对幅度: %.3f %.3f %.3f %.3f %.3f\n', ...
    sorted(1:5)/sorted(1));
end
