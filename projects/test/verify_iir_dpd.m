function verify_iir_dpd()
%VERIFY_IIR_DPD 用 FIR 逼近 IIR 滤波器系数并验证其能否支撑 DPD（端到端可行性测试）
%
% 流程（对应 7 个步骤）：
%   1 导入数据        x / y0（输入 / 无 DPD 的 PA 输出）
%   2 生成 IIR 椭圆低通滤波器（带宽 80 MHz）
%   3 滤波            y0 -> y_f（带限观测，模拟采集路径，循环卷积）
%   4 提取 FIR 系数   用 FIR 逼近 IIR，产出滤波器系数 h
%   5 用 BLMP 训练 PA 模型
%   6 估计全带输出
%   7 做 DPD 并评估 NMSE / EVM

%% 1. 导入数据（x: 基带输入, y0: PA 输出, 单位见数据集说明）
d = load(fullfile('datasets', '测试', '2026.6.10', 'test1', 'x_v0.6_v1.mat'));
x = struct2array(d); x = x(:);
d = load(fullfile('datasets', '测试', '2026.6.10', 'test1', 'y0_v0.6_v1.mat'));
y0 = struct2array(d); y0 = y0(:);

%% 2. 生成 IIR 椭圆低通滤波器（带宽 80 MHz）
Fs = 1000;          % 采样率 MHz
BW = 50;            % 带宽 MHz（单边截止）
Rp = 0.5;           % 通带纹波 dB
Rs = 80;            % 阻带衰减 dB
Wp = BW / (Fs/2);               % 通带边缘（归一化）
Ws = (BW + 3) / (Fs/2);         % 阻带边缘（归一化）
[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn, 'low');
fprintf('IIR 椭圆低通：阶数 n=%d, 带宽 BW=%d MHz\n', n, BW);

%% 3. 滤波：y0 对齐到 x 后过椭圆滤波器，得到带限观测 y_f（循环卷积）
%[y0a, ~, ~] = align_circular_signal(y0, x);
y0a = y0;
y_f = circfilter(y0a, b, a);

%% 4. 提取 FIR 系数（用 FIR 逼近 IIR，基于 F_matrix 循环卷积矩阵）
% 最终产出 FIR 滤波器系数 h（脉冲响应），用有限阶 FIR 逼近 IIR。
% 扫描 FIR 长度 L，按重构 NMSE 选最优。
L_list = [64 128 256 512 1024 2048];
nmse_fit = zeros(size(L_list));
for i = 1:numel(L_list)
    h_i = extract_fir_coeff(y0a, y_f, L_list(i));
    nmse_fit(i) = get_NMSE(circfilter(y0a, h_i, 1), y_f);
    fprintf('L=%3d: 重构 NMSE = %.2f dB\n', L_list(i), nmse_fit(i));
end
[~, i_best] = min(nmse_fit);
L = L_list(i_best);
h = extract_fir_coeff(y0a, y_f, L);
y_f_hat = circfilter(y0a, h, 1);
fprintf('选用 FIR 长度 L=%d（IIR 阶数 n=%d）\n', L, n);

%% 4b. 提取验证对比图：y / y_f / y_f_hat 功率谱
nfft = 4096;
win = hann(nfft);
[p0, f] = pwelch(y0a, win, nfft/2, nfft, Fs);
[p1, ~] = pwelch(y_f, win, nfft/2, nfft, Fs);
[p2, ~] = pwelch(y_f_hat, win, nfft/2, nfft, Fs);
figure;
plot(f - Fs/2, 10*log10(fftshift(p0)), 'k-', 'LineWidth', 1.2); hold on;
plot(f - Fs/2, 10*log10(fftshift(p1)), 'b-', 'LineWidth', 1.2);
plot(f - Fs/2, 10*log10(fftshift(p2)), 'g-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('功率谱密度 (dB)');
legend('y (滤波器输入 y0a)', 'y_f (IIR 真值)', 'y_f_{hat} (FIR 逼近)');
title('FIR 逼近 IIR 提取验证');
grid on;
%% 5. 用 BLMP 训练 PA 模型
% BLMP = 全带 MP 基函数过反馈 FIR 滤波器（此处即提取的 FIR h）。
% circfilter(X,h,1) 等价于 geth_even_BLMP(x,K,M,h) 的循环卷积版。
K = 7; M = 5;                   % 非线性阶数 / 记忆深度
y_full = get_MG(K, M, 3, x, h, y_f); 
%% 6. 估计全带输出
fprintf('全带估计 NMSE(y_full, y0): %.2f dB\n', get_NMSE(y_full, y0));
%% 6b. 全带估计对比图：y0a / y_f / y_full 功率谱
[p3, f] = pwelch(y0a, win, nfft/2, nfft, Fs);
[p4, ~] = pwelch(y_f, win, nfft/2, nfft, Fs);
[p5, ~] = pwelch(y_full, win, nfft/2, nfft, Fs);
figure;
plot(f - Fs/2, 10*log10(fftshift(p3)), 'k-', 'LineWidth', 1.2); hold on;
plot(f - Fs/2, 10*log10(fftshift(p4)), 'b-', 'LineWidth', 1.2);
plot(f - Fs/2, 10*log10(fftshift(p5)), 'g-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('功率谱密度 (dB)');
legend('y0a (真实 PA 输出)', 'y_f (带限观测)', 'y_{full} (估计全带)');
title('全带输出估计验证');
grid on;
end

function h = extract_fir_coeff(x, y_f, L)
%EXTRACT_FIR_COEFF 基于 F_matrix 的循环卷积矩阵提取 FIR 系数（脉冲响应）。
% 用有限阶 FIR 逼近 IIR：y_f ≈ 循环卷积(x, h)，h 为 L 阶 FIR 系数。
% C(·) 由 F_matrix 的 F + [S;zeros] 得到（线性卷积 + 循环边界修正）。
    N = numel(x);
    [Fx, Sx] = F_matrix(x, L);
    Cx = Fx + [Sx; zeros(N - (L - 1), L)];
    % 带限信号下卷积矩阵病态，加小 ridge 正则化稳定数值，避免 NaN
    lambda = 1e-6 * norm(Cx, 'fro')^2 / size(Cx, 2);
    h = (Cx' * Cx + lambda * eye(L)) \ (Cx' * y_f);
end

function y = circfilter(x, b, a)
%CIRCFILTER 循环卷积滤波：Y(f) = H(f)·X(f)，无瞬态（循环信号 y(-1)=y(N)）。
% 对矩阵按列滤波；FIR 时 a=1。
    N = size(x, 1);
    H = freqz(b, a, N, 'whole');
    y = ifft(fft(x) .* H);
end


