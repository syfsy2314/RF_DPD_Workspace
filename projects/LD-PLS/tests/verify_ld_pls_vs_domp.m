function results = verify_ld_pls_vs_domp()
%VERIFY_LD_PLS_VS_DOMP 比较全带宽能量惩罚 LD-PLS 与 DOMP。
% 两种方法使用相同的 H1+H2 测量链路：前 60% 拟合，中间 20% 给 DOMP
% 选阶，末尾 20% 评估。LD-PLS 固定 mu=1，不扫描惩罚参数。

%% 参数和数据
projectDir = fileparts(fileparts(mfilename('fullpath')));
dataX = load(fullfile(projectDir, 'data', 'x.mat'), 'x');
dataY = load(fullfile(projectDir, 'data', 'y0.mat'), 'y0');
x = dataX.x(:);
y0 = dataY.y0(:);

fs = 1000;  % 采样率，MHz；下面的带宽和频率偏移也以 MHz 为单位。
bandwidths = [70 60 50 40];
K = 7;
M = 5;
lambda = 10;
mu = 1;  % 固定全带宽能量惩罚强度，不再逐带宽扫描。
T_delta = 0;  % 关闭边际增益条件；仅由得分矩阵的条件数决定停止。
kappa_G_max = 10^4.5;
maxO = round(0.5 * K * (M + 1));  % 最多提取 21 个潜变量。

nTrain = floor(0.6 * numel(x));
nValidation = floor(0.2 * numel(x));
train = 1:nTrain;
validation = nTrain + (1:nValidation);
test = (nTrain + nValidation + 1):numel(x);

%% 全带宽基函数
X = geth_even_1d(x, K, M);
nBandwidths = numel(bandwidths);
ldO = zeros(nBandwidths, 1);
ldValidationNmse = zeros(nBandwidths, 1);
ldTestNmse = zeros(nBandwidths, 1);
dompOrder = zeros(nBandwidths, 1);
dompValidationNmse = zeros(nBandwidths, 1);
dompTestNmse = zeros(nBandwidths, 1);

figure('Name', '全带宽能量惩罚 LD-PLS 与 DOMP');

for ib = 1:nBandwidths
    BW = bandwidths(ib);

    %% H1 + H2：模拟实际链路
    h1 = design_filter([BW + 5 BW + 10], 0.5, 1e-4, fs);
    h2 = design_filter([BW BW + 3], 0.1, 1e-6, fs);
    y_f = filter(h1, 1, y0);
    [y_f, ~, ~] = align_circular_signal(y_f, x);
    y_bf = filter(h2, 1, y_f);
    % H_x 是模型通过 H2 后的观测基函数；X 用于重建全带宽输出。
    H_x = filter(h2, 1, X);

    %% LD-PLS：固定 mu=1，自适应选择潜变量数
    [hLD, ldO(ib)] = get_energy_regularized_LD_PLS( ...
        X(train, :), H_x(train, :), y_bf(train), lambda, ...
        mu, T_delta, kappa_G_max, maxO);
    ldValidationNmse(ib) = get_NMSE(X(validation, :) * hLD, y0(validation));
    yLD = X(test, :) * hLD;
    ldTestNmse(ib) = get_NMSE(yLD, y0(test));

    %% DOMP：在相同训练段拟合，在相同验证段选阶
    labels = cell(1, size(H_x, 2));
    % hHistory 的每一列对应一个 DOMP 阶数的系数，可直接作用于 X。
    [~, ~, ~, hHistory] = omp_domp( ...
        H_x(train, :), y_bf(train), labels, 'DOMP');
    validationByOrder = zeros(1, size(hHistory, 2));
    for order = 1:size(hHistory, 2)
        yValidation = X(validation, :) * hHistory(:, order);
        validationByOrder(order) = get_NMSE(yValidation, y0(validation));
    end
    [dompValidationNmse(ib), dompOrder(ib)] = min(validationByOrder);
    yDOMP = X(test, :) * hHistory(:, dompOrder(ib));
    dompTestNmse(ib) = get_NMSE(yDOMP, y0(test));

    %% 测试段功率谱：与 main.m 一样按各信号均方根归一化
    subplot(nBandwidths, 1, ib);
    plot_psd(y0(test), fs, 'k-');
    hold on;
    plot_psd(yLD, fs, 'm-');
    plot_psd(yDOMP, fs, 'b-');
    hold off;
    xlim([-120 120]);
    grid on;
    title(sprintf('%d MHz: LD-PLS %.2f dB (O=%d, mu=%g), DOMP %.2f dB (order=%d)', ...
        BW, ldTestNmse(ib), ldO(ib), mu, dompTestNmse(ib), dompOrder(ib)));
    xlabel('频率偏移 (MHz)');
    ylabel('功率谱密度 (dB)');
    if ib == 1
        legend('PA 输出', 'LD-PLS', 'DOMP', 'Location', 'best');
    end
end

%% 数值结果
fprintf('H1+H2；训练/验证/测试 = 60/20/20%%；mu=%g，T_delta=%g，kappa_G_max=%g\n', ...
    mu, T_delta, kappa_G_max);
fprintf('BW(MHz)  LD: O / mu / 验证 / 测试(dB)       DOMP: 阶数 / 验证 / 测试(dB)\n');
for ib = 1:nBandwidths
    fprintf('%3d       %2d / %-5g / %7.2f / %7.2f          %2d / %7.2f / %7.2f\n', ...
        bandwidths(ib), ldO(ib), mu, ldValidationNmse(ib), ...
        ldTestNmse(ib), dompOrder(ib), dompValidationNmse(ib), ...
        dompTestNmse(ib));
end

results.bandwidths = bandwidths;
results.ldO = ldO;
results.ldMu = mu;
results.ldValidationNmse = ldValidationNmse;
results.ldTestNmse = ldTestNmse;
results.dompOrder = dompOrder;
results.dompValidationNmse = dompValidationNmse;
results.dompTestNmse = dompTestNmse;
results.T_delta = T_delta;
results.kappa_G_max = kappa_G_max;
end

function h = design_filter(edges, rippleDb, stopRipple, fs)
%DESIGN_FILTER 按通带/阻带边缘设计 FIR 低通滤波器。
% edges=[通带边缘, 阻带边缘]，单位 MHz；fs 单位 MHz。
% rippleDb 为通带波纹 dB，stopRipple 为线性幅度阻带容差；h 为 FIR 系数。
passRipple = (10^(rippleDb / 20) - 1) / (10^(rippleDb / 20) + 1);
[order, frequency, amplitude, weights] = firpmord( ...
    edges / (fs / 2), [1 0], [passRipple stopRipple]);
h = firpm(order, frequency, amplitude, weights, {20});
end

function plot_psd(y, fs, lineStyle)
%PLOT_PSD 绘制均方根归一化后的双边功率谱，横轴单位 MHz。
% 乘以 100 沿用 main.m 的纵轴缩放；不同曲线比较相对谱形和泄漏。
nfft = 1024;
window = hann(nfft);
[power, frequency] = pwelch(y / rms(y), window, nfft / 2, nfft, fs);
plot(frequency - fs / 2, 10 * log10(fftshift(power * 100)), lineStyle);
end
