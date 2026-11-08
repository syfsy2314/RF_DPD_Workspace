function diag_phase_coverage()
% 量化：LD-PLS 偶数阶 BLMP 基函数（x·|x|^(k-1)，纯同相）对相位（正交）失真的覆盖能力。
% 结论应是：能覆盖 AM/AM，但覆盖不了 AM/PM 与 OMUX 相位。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d);
x = x(:);
K = 7; M = 5;

% 幅度失真系数（偶数阶，可被基函数精确表达）
am_coef = 1 + 0.3*abs(x).^2 + 0.1*abs(x).^4;

% 1) 纯幅度失真（AM/AM only，相位 = angle(x)）
y_am = x .* am_coef;

% 2) 幅度 + 相位失真（AM/AM + AM/PM，相位随幅度变化）
y_ampm = y_am .* exp(1j*0.5*abs(x).^2);

% 3) 幅度 + 频率相关相位（AM/AM + OMUX 全通相位）
r = 0.8; theta = 0.6*pi;
b = [r^2, -2*r*cos(theta), 1]; a = [1, -2*r*cos(theta), r^2];
y_omux = filter(b, a, y_am);

% 偶数阶 BLMP 基函数（同相）
H = geth_even_1d(x, K, M);

% 复数最小二乘拟合
y_hat_am    = H * ((H'*H) \ (H'*y_am));
y_hat_ampm  = H * ((H'*H) \ (H'*y_ampm));
y_hat_omux  = H * ((H'*H) \ (H'*y_omux));

fprintf('--- 偶数阶 BLMP 基函数覆盖能力 ---\n');
fprintf('纯幅度失真 AM/AM         : NMSE = %6.2f dB\n', get_NMSE(y_hat_am,   y_am));
fprintf('幅度+AM/PM 相位失真      : NMSE = %6.2f dB\n', get_NMSE(y_hat_ampm, y_ampm));
fprintf('幅度+OMUX 频率相位失真   : NMSE = %6.2f dB\n', get_NMSE(y_hat_omux, y_omux));
end
