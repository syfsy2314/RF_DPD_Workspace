function test_phase_budget()
%TEST_PHASE_BUDGET 相位预算观测实验：三组相位差随频率的变化
%  ① x vs y   ：PA 造成
%  ② y vs y_f ：滤波器造成（椭圆 IIR，非线性相位）
%  ③ x vs y_f ：PA + 滤波器共同造成
%  自洽校验：③ = ① + ②（恒等）；② = angle(H1)（理论）。

%% 参数
Fs = 1000;          % 采样率 MHz
Rp = 0.5;           % H1 通带纹波 dB
Rs = 80;            % H1 阻带衰减 dB

%% 数据
d = load(fullfile('datasets', '测试', '2026.6.10', 'test1', 'x_v0.6_v1.mat'));
x = struct2array(d); x = x(:);
d = load(fullfile('datasets', '测试', '2026.6.10', 'test1', 'y0_v0.6_v1.mat'));
y0 = struct2array(d); y0 = y0(:);
N = numel(x);

%% 估计信号带边（99% 能量），让 H1 过渡带落在信号带边
f = fftshift((-N/2:N/2-1).') / N;       % 归一化频率
fMHz = f * Fs;                          % MHz
X = fftshift(fft(x));
cdf = cumsum(abs(X).^2) / sum(abs(X).^2);
i_lo = find(cdf >= 0.005, 1);
i_hi = find(cdf >= 0.995, 1);
f_edge = max(abs(f(i_lo)), abs(f(i_hi)));
fprintf('信号带边 = ±%.4f (%.1f MHz)\n', f_edge, f_edge * Fs);

%% H1：椭圆低通 IIR（非线性相位）
Wp = 0.8 * f_edge;
Ws = 1.5 * f_edge;
[N1, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b1, a1] = ellip(N1, Rp, Rs, Wn, 'low');
fprintf('H1(椭圆 IIR) 阶数 = %d, 通带 %.4f / 阻带 %.4f\n', N1, Wp, Ws);

%% 对齐 + 链路：y0 对齐到 x，再经 H1 得到带限输出 y_f
y0a = align_to(x, y0);
y_f = filter(b1, a1, y0a);

%% 频域相位差（主值，rad）
Y0 = fftshift(fft(y0a));
Yf = fftshift(fft(y_f));
phi_xy = wrapToPi(angle(Y0) - angle(X));   % ① PA
phi_yf = wrapToPi(angle(Yf) - angle(Y0));  % ② 滤波器
phi_xf = wrapToPi(angle(Yf) - angle(X));   % ③ 总

%% 掩码：x 能量高于峰值 -40 dB 的频段（带外相位差是噪声）
mask = abs(X).^2 > max(abs(X).^2) * 1e-4;

%% 理论对照 + 加性校验
H1f = fftshift(freqz(b1, a1, N, 'whole'));
phi_H1 = angle(H1f);
err_theory = wrapToPi(phi_yf - phi_H1);            % ② vs angle(H1)
err_add    = wrapToPi(phi_xf - (phi_xy + phi_yf)); % ③ vs ①+②
fprintf('--- 校验（带内 RMS，rad）---\n');
fprintf('② 滤波器 vs angle(H1): %.3e rad\n', rms(err_theory(mask)));
fprintf('③ 总相位 vs ①+②    : %.3e rad\n', rms(err_add(mask)));

%% 对比图：三组相位差随频率的变化
figure(1); clf;
fl = [-f_edge f_edge] * Fs;

subplot(3,1,1);
plot(fMHz(mask), phi_xy(mask)*180/pi, 'b.', 'MarkerSize', 3);
xlabel('频率 (MHz)'); ylabel('相位差 (deg)');
title('① x 与 y：PA 造成的相位差'); grid on; xlim(fl);

subplot(3,1,2);
plot(fMHz(mask), phi_yf(mask)*180/pi, 'r.', 'MarkerSize', 3); hold on;
plot(fMHz(mask), phi_H1(mask)*180/pi, 'k-', 'LineWidth', 1);
xlabel('频率 (MHz)'); ylabel('相位差 (deg)');
title('② y 与 y_f：滤波器造成的相位差'); grid on;
legend({'观测', 'angle(H1) 理论'}, 'Location', 'best'); xlim(fl);

subplot(3,1,3);
plot(fMHz(mask), phi_xf(mask)*180/pi, 'k.', 'MarkerSize', 3); hold on;
plot(fMHz(mask), wrapToPi(phi_xy(mask) + phi_yf(mask))*180/pi, ...
    'm.', 'MarkerSize', 3);
xlabel('频率 (MHz)'); ylabel('相位差 (deg)');
title('③ x 与 y_f：PA + 滤波器共同造成的相位差'); grid on;
legend({'观测', '①+② 加性验证'}, 'Location', 'best'); xlim(fl);
end

%% 局部函数
function ya = align_to(x, y)
% 整数时延 + 常数相位对齐（参考 LD-PLS verify_routeB）
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
