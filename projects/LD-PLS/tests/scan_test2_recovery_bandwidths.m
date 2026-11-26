function results = scan_test2_recovery_bandwidths(inputIqPath, bandwidthsMHz)
%SCAN_TEST2_RECOVERY_BANDWIDTHS 对比 test2 的两种 LD-PLS 和 DOMP。
% 从仓库根目录执行 startup，再传入配对的原始交错 int16 IQ 文件。
% bandwidthsMHz 是 H2 的单边通带边缘，默认扫描 30 至 80 MHz。
% 用同一组全带 y0 扫描选阶并评价，结果为样本内全带 NMSE。

if nargin < 2
    bandwidthsMHz = 30:10:80;
end

fsMHz = 1000;
validateattributes(bandwidthsMHz, {'numeric'}, ...
    {'vector', 'nonempty', 'real', 'finite', 'positive', '<', fsMHz / 2 - 10});
bandwidthsMHz = bandwidthsMHz(:);

repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dataDir = fullfile(repoRoot, 'datasets', '测试', '2026.11.24', 'test2');
data = load(fullfile(dataDir, 'y0.mat'), 'y0');
y0 = data.y0(:);
x = read_input_iq(inputIqPath, numel(y0));
x = x / max(abs(x));

K = 7;
M = 5;
lambda = 10;
X = geth_even_1d(x, K, M);
nFeatures = size(X, 2);
muValues = [1, 0];
ldBestNmse = zeros(numel(bandwidthsMHz), numel(muValues));
ldBestOrder = zeros(size(ldBestNmse));
dompBestNmse = zeros(size(bandwidthsMHz));
dompBestOrder = zeros(size(bandwidthsMHz));

for idx = 1:numel(bandwidthsMHz)
    bw = bandwidthsMHz(idx);
    % 与 LD-PLS main.m 使用相同的 FIR 设计参数和对齐顺序。
    h1PassRipple = (10^(0.5 / 20) - 1) / (10^(0.5 / 20) + 1);
    h1 = design_filter([bw + 5, bw + 10], h1PassRipple, 1e-4, fsMHz);
    h2 = design_filter([bw, bw + 3], 0.0057501127785, 1e-6, fsMHz);
    yH1 = align_circular_signal(filter(h1, 1, y0), x);
    yBand = filter(h2, 1, yH1);
    XBand = filter(h2, 1, X);

    for im = 1:numel(muValues)
        scores = zeros(1, nFeatures - 1);
        for order = 1:numel(scores)
            a = get_energy_regularized_LD_PLS( ...
                X, XBand, yBand, lambda, muValues(im), 0, inf, order);
            scores(order) = get_NMSE(X * a, y0);
        end
        [ldBestNmse(idx, im), ldBestOrder(idx, im)] = min(scores);
    end

    [~, ~, ~, dompHistory] = omp_domp(XBand, yBand, cell(1, nFeatures), 'DOMP');
    scores = zeros(1, size(dompHistory, 2));
    for order = 1:numel(scores)
        scores(order) = get_NMSE(X * dompHistory(:, order), y0);
    end
    [dompBestNmse(idx), dompBestOrder(idx)] = min(scores);

    fprintf('BW %g MHz | LD-PLS mu=1 %.3f dB (O=%d) | mu=0 %.3f dB (O=%d) | DOMP %.3f dB (%d terms)\n', ...
        bw, ldBestNmse(idx, 1), ldBestOrder(idx, 1), ...
        ldBestNmse(idx, 2), ldBestOrder(idx, 2), dompBestNmse(idx), dompBestOrder(idx));
end

results = table(bandwidthsMHz, ldBestNmse(:, 1), ldBestOrder(:, 1), ...
    ldBestNmse(:, 2), ldBestOrder(:, 2), dompBestNmse, dompBestOrder, ...
    'VariableNames', {'BW_MHz', 'LDPLS_Mu1_NMSE_dB', 'LDPLS_Mu1_Order', ...
    'LDPLS_Mu0_NMSE_dB', 'LDPLS_Mu0_Order', 'DOMP_NMSE_dB', 'DOMP_Terms'});
end

function h = design_filter(edgesMHz, passRipple, stopRipple, fsMHz)
[order, frequencies, amplitudes, weights] = firpmord( ...
    edgesMHz / (fsMHz / 2), [1 0], [passRipple, stopRipple]);
h = firpm(order, frequencies, amplitudes, weights, {20});
end

function x = read_input_iq(inputIqPath, expectedSamples)
fid = fopen(inputIqPath, 'rb', 'ieee-le');
if fid < 0
    error('scan_test2_recovery_bandwidths:InputFile', 'Cannot open IQ file: %s', inputIqPath);
end
closeFile = onCleanup(@() fclose(fid));
[raw, count] = fread(fid, [2, inf], 'int16=>double');
if count ~= 2 * expectedSamples
    error('scan_test2_recovery_bandwidths:SampleCount', ...
        'Expected %d IQ samples to match test2 y0, but found %g.', expectedSamples, count / 2);
end
x = complex(raw(1, :), raw(2, :)).';
end
