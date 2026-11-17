function results = verify_ld_pls_vs_domp()
%VERIFY_LD_PLS_VS_DOMP 比较 DOMP、原 LD-PLS 和全带惩罚 LD-PLS。
% 五组方法使用项目初始数据及相同 H1+H2 链路，全部样本用于拟合
% 和全带 NMSE 评价。DOMP 扫描选阶；两种 LD-PLS 均比较自适应与扫描 O。
% 拟合、选阶和评价共用样本，仅反映样本内效果。

%% 参数和数据
projectDir = fileparts(fileparts(mfilename('fullpath')));
dataDir = fullfile(projectDir, 'data');
dataX = load(fullfile(dataDir, 'x.mat'), 'x');
dataY = load(fullfile(dataDir, 'y0.mat'), 'y0');
x = dataX.x(:);
y0 = dataY.y0(:);

fs = 1000;  % 采样率，MHz；下面的带宽和频率偏移也以 MHz 为单位。
bandwidths = [80 70 60 50 40 30];  % 修改此列表即可调整比较带宽。
K = 7;
M = 5;
lambda = 10;
mu = 0.1;  % 六组带宽共用同一惩罚强度。
T_delta = 0;  % 关闭边际增益条件；仅由得分矩阵的条件数决定停止。
kappa_G_max = 10^4.5;
maxO = round(0.5 * K * (M + 1));  % 最多提取 21 个潜变量。

%% 全带宽基函数
X = geth_even_1d(x, K, M);
nBandwidths = numel(bandwidths);
ldO = zeros(nBandwidths, 1);
ldNmse = zeros(nBandwidths, 1);
ldScanO = zeros(nBandwidths, 1);
ldScanNmse = zeros(nBandwidths, 1);
ldNmseByO = zeros(nBandwidths, maxO);
energyO = zeros(nBandwidths, 1);
energyNmse = zeros(nBandwidths, 1);
energyScanO = zeros(nBandwidths, 1);
energyScanNmse = zeros(nBandwidths, 1);
energyNmseByO = zeros(nBandwidths, maxO);
dompOrder = zeros(nBandwidths, 1);
dompNmse = zeros(nBandwidths, 1);

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

    %% 原 LD-PLS：mu=0，潜变量数由停止准则决定
    [hLD, ldO(ib)] = get_energy_regularized_LD_PLS( ...
        X, H_x, y_bf, lambda, 0, T_delta, kappa_G_max, maxO);
    yLD = X * hLD;
    ldNmse(ib) = get_NMSE(yLD, y0);

    %% 原 LD-PLS：关闭自适应停止，扫描固定 O=1:maxO
    for candidateO = 1:maxO
        [hScanned, fittedO] = get_energy_regularized_LD_PLS( ...
            X, H_x, y_bf, lambda, 0, 0, inf, candidateO);
        assert(fittedO == candidateO, '原 LD-PLS 扫描 O 时提前停止');
        ldNmseByO(ib, candidateO) = get_NMSE(X * hScanned, y0);
    end
    [ldScanNmse(ib), ldScanO(ib)] = min(ldNmseByO(ib, :));

    %% 全带能量惩罚 LD-PLS：固定 mu=0.1，自适应选潜变量
    [hEnergy, energyO(ib)] = get_energy_regularized_LD_PLS( ...
        X, H_x, y_bf, lambda, mu, T_delta, kappa_G_max, maxO);
    yEnergy = X * hEnergy;
    energyNmse(ib) = get_NMSE(yEnergy, y0);

    %% 全带能量惩罚 LD-PLS：关闭自适应停止，扫描固定 O=1:maxO
    % 每个 O 都以相同全信号的全带 NMSE 评分，作为样本内选阶上限。
    for candidateO = 1:maxO
        [hScanned, fittedO] = get_energy_regularized_LD_PLS( ...
            X, H_x, y_bf, lambda, mu, 0, inf, candidateO);
        assert(fittedO == candidateO, '扫描 O 时提前停止');
        energyNmseByO(ib, candidateO) = get_NMSE(X * hScanned, y0);
    end
    [energyScanNmse(ib), energyScanO(ib)] = min(energyNmseByO(ib, :));

    %% DOMP：全信号拟合，并用同一全信号的全带 NMSE 扫描选阶
    labels = cell(1, size(H_x, 2));
    % hHistory 的每一列对应一个 DOMP 阶数的系数，可直接作用于 X。
    [~, ~, ~, hHistory] = omp_domp(H_x, y_bf, labels, 'DOMP');
    nmseByOrder = zeros(1, size(hHistory, 2));
    for order = 1:size(hHistory, 2)
        nmseByOrder(order) = get_NMSE(X * hHistory(:, order), y0);
    end
    [dompNmse(ib), dompOrder(ib)] = min(nmseByOrder);
    fprintf('%d MHz 完成：LD-PLS 自适应 O=%d / %.3f dB；扫描 O=%d / %.3f dB；惩罚版自适应 O=%d / %.3f dB；扫描 O=%d / %.3f dB\n', ...
        BW, ldO(ib), ldNmse(ib), ldScanO(ib), ldScanNmse(ib), ...
        energyO(ib), energyNmse(ib), energyScanO(ib), energyScanNmse(ib));
end

%% 全信号全带恢复 NMSE：数值表与带宽曲线
fprintf('项目初始数据；H1+H2；全部信号训练并评价（样本内）；mu=%g，lambda=%g，T_delta=%g，kappa_G_max=%g\n', ...
    mu, lambda, T_delta, kappa_G_max);
resultsTable = table(bandwidths(:), dompOrder, ldO, ldScanO, energyO, energyScanO, ...
    dompNmse, ldNmse, ldScanNmse, energyNmse, energyScanNmse, ...
    'VariableNames', {'BandwidthMHz', 'DOMP_Order', 'LD_Adaptive_O', ...
    'LD_Scanned_O', 'Energy_Adaptive_O', 'Energy_Scanned_O', ...
    'DOMP_NMSE_dB', 'LD_Adaptive_NMSE_dB', 'LD_Scanned_NMSE_dB', ...
    'Energy_Adaptive_NMSE_dB', 'Energy_Scanned_NMSE_dB'});
disp(resultsTable);

[sortedBandwidths, order] = sort(bandwidths);
figure('Name', '五组方法的全带恢复 NMSE');
plot(sortedBandwidths, dompNmse(order), 'b-o', ...
    sortedBandwidths, ldNmse(order), 'g-s', ...
    sortedBandwidths, ldScanNmse(order), 'c-v', ...
    sortedBandwidths, energyNmse(order), 'm-^', ...
    sortedBandwidths, energyScanNmse(order), 'r-d', 'LineWidth', 1.5);
grid on;
xlabel('观测滤波器带宽 (MHz)');
ylabel('全信号全带恢复 NMSE (dB)');
legend('DOMP（扫描选阶）', 'LD-PLS（自适应 O）', 'LD-PLS（扫描 O）', ...
    '全带惩罚 LD-PLS（自适应 O）', '全带惩罚 LD-PLS（扫描 O）', ...
    'Location', 'best');

results.bandwidths = bandwidths;
results.dataDir = dataDir;
results.ldO = ldO;
results.ldNmse = ldNmse;
results.ldScanO = ldScanO;
results.ldScanNmse = ldScanNmse;
results.ldNmseByO = ldNmseByO;
results.energyO = energyO;
results.energyMu = mu;
results.energyNmse = energyNmse;
results.energyScanO = energyScanO;
results.energyScanNmse = energyScanNmse;
results.energyNmseByO = energyNmseByO;
results.dompOrder = dompOrder;
results.dompNmse = dompNmse;
results.T_delta = T_delta;
results.kappa_G_max = kappa_G_max;
results.table = resultsTable;
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
