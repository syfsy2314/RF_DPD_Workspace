%% 开关 2 同时训练相位和幅度，再切回开关 1 恢复 PA 全带输出
% 延续 Phase_Calibration 的信号、互谱、LD-PLS 与评分设置。
% 开关 2 仅用原始输入 x 和旁路反馈 x_f；H1 的通带波纹改为 1 dB。
clear;
close all;

projectDir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(projectDir));
dataDir = fullfile(root, 'datasets', '测试', '2026.6.10', 'test1');
resultDir = fullfile(projectDir, 'results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

Fs = 1000;                     % MHz
RpDb = 1;                      % H1 通带波纹，dB
RsDb = 80;
H1passMHz = 85;
H1stopMHz = 90;
H2edgesMHz = [80, 44];
Dpass2 = 0.0057501127785;
Dstop2 = 1e-4;
K = 7;
M = 5;
S = 0.3;
lamda = 1;
guard = 5000;
segmentLength = 4096;
overlap = segmentLength / 2;
nfft = 8192;
smoothingBins = 5;
phaseTapCount = 4001;
taperMHz = 5;

d = load(fullfile(dataDir, 'x_v0.6_v1.mat'));
x = d.x(:);
d = load(fullfile(dataDir, 'y0_v0.6_v1.mat'));
y0 = d.y0(:);

%% 同一条 H1 作用于开关 1 的 PA 输出和开关 2 的原始输入
[h1Order, wn] = ellipord(H1passMHz / (Fs / 2), ...
    H1stopMHz / (Fs / 2), RpDb, RsDb);
[b1, a1] = ellip(h1Order, RpDb, RsDb, wn, 'low');
yFilteredRaw = filter(b1, a1, y0);
[yFiltered, xAligned] = align_circular_signal(yFilteredRaw, x);
xFilteredRaw = filter(b1, a1, x);
xFiltered = align_circular_signal(xFilteredRaw, xAligned);
X = geth_even_1d(xAligned, K, M);

%% 相位校准原实验的互谱设置与可靠频点选择
idx = (guard + 1):(numel(xAligned) - guard);
window = hann(segmentLength, 'periodic');
options = {window, overlap, nfft, Fs, 'centered'};
Syx = movmean(cpsd(yFiltered(idx), xAligned(idx), options{:}), ...
    smoothingBins);
Sxx = movmean(pwelch(xAligned(idx), options{:}), smoothingBins);
Syy = movmean(pwelch(yFiltered(idx), options{:}), smoothingBins);
coherence = min(1, abs(Syx).^2 ./ (Sxx .* Syy));
fGrid = (-nfft / 2:nfft / 2 - 1).' * Fs / nfft;
reliable = abs(fGrid) <= 80 & ...
    10 * log10(Sxx / max(Sxx)) >= -30 & coherence >= 0.8;
fReliable = fGrid(reliable);
omegaReliable = 2 * pi * fReliable / Fs;
design = [ones(nnz(reliable), 1), omegaReliable];

% 开关 2 的 x_f=H1(x) 同时提供相位与幅度观测。
Sxfx = movmean(cpsd(xFiltered(idx), xAligned(idx), options{:}), ...
    smoothingBins);
transfer = Sxfx(reliable) ./ Sxx(reliable);
phaseObserved = unwrap(angle(transfer));
phaseResidual = phaseObserved - design * (design \ phaseObserved);
magnitudeObservedDb = 20 * log10(abs(transfer));
[~, centerIndex] = min(abs(fReliable));
magnitudeDb = magnitudeObservedDb - magnitudeObservedDb(centerIndex);

% 两组只在模型侧算子上不同；真实反馈与评分样本完全相同。
hPhase = create_calibration_fir(fReliable, phaseResidual, ...
    zeros(size(magnitudeDb)), Fs, nfft, phaseTapCount, taperMHz);
hPhaseAmplitude = create_calibration_fir(fReliable, phaseResidual, ...
    magnitudeDb, Fs, nfft, phaseTapCount, taperMHz);

%% 离线核对：真实 H1 只用于生成观测、参考和评分
h1Impulse = filter(b1, a1, [1; zeros(phaseTapCount - 1, 1)]);
trueResponse = freqz(b1, a1, fReliable, Fs);
trueMagnitudeDb = 20 * log10(abs(trueResponse));
trueMagnitudeDb = trueMagnitudeDb - trueMagnitudeDb(centerIndex);
truePhase = unwrap(angle(trueResponse));
truePhaseResidual = truePhase - design * (design \ truePhase);
magnitudeErrorRmsDb = rms(magnitudeDb - trueMagnitudeDb);
phaseErrorRmsRad = rms(phaseResidual - truePhaseResidual);

methodNames = {'仅 H2', '开关 2 相位', '开关 2 相位+幅度', '已知 H1 参考'};
nmseDb = zeros(numel(H2edgesMHz), numel(methodNames));
h2Orders = zeros(numel(H2edgesMHz), 1);
for iBw = 1:numel(H2edgesMHz)
    BW = H2edgesMHz(iBw);
    [h2Orders(iBw), fo, ao, weights] = firpmord( ...
        [BW, BW + 3] / (Fs / 2), [1 0], [Dpass2, Dstop2]);
    h2 = firpm(h2Orders(iBw), fo, ao, weights, {20});
    yFeedback = filter(h2, 1, yFiltered);
    filters = {h2(:), conv(hPhase, h2(:)), ...
        conv(hPhaseAmplitude, h2(:)), conv(h1Impulse, h2(:))};
    for iMethod = 1:numel(methodNames)
        hModel = filters{iMethod};
        reference = filter(hModel, 1, xAligned);
        feedbackAligned = align_circular_signal(yFeedback, reference);
        yEstimate = get_LD_PLS(K, M, X, xAligned, hModel, ...
            feedbackAligned, S, lamda);
        yEstimate = align_circular_signal(yEstimate, y0);
        nmseDb(iBw, iMethod) = get_NMSE(yEstimate(idx), y0(idx));
        if BW == 44 && iMethod == 2
            yPhase44 = yEstimate;
        elseif BW == 44 && iMethod == 3
            yPhaseAmplitude44 = yEstimate;
        end
    end
end

%% H1 全通带响应：观测曲线只保留可靠频点，模型曲线显示带外延拓
h1FrequencyMHz = fGrid(abs(fGrid) <= H1passMHz);
h1Reliable = ismember(h1FrequencyMHz, fReliable);
h1Design = [ones(numel(h1FrequencyMHz), 1), ...
    2 * pi * h1FrequencyMHz / Fs];
h1TrueResponse = freqz(b1, a1, h1FrequencyMHz, Fs);
h1ModelResponse = freqz(hPhaseAmplitude, 1, h1FrequencyMHz, Fs);
[~, h1CenterIndex] = min(abs(h1FrequencyMHz));
h1TrueMagnitudeDb = 20 * log10(abs(h1TrueResponse));
h1TrueMagnitudeDb = h1TrueMagnitudeDb - ...
    h1TrueMagnitudeDb(h1CenterIndex);
h1ModelMagnitudeDb = 20 * log10(abs(h1ModelResponse));
h1ModelMagnitudeDb = h1ModelMagnitudeDb - ...
    h1ModelMagnitudeDb(h1CenterIndex);
h1TruePhase = unwrap(angle(h1TrueResponse));
h1ModelPhase = unwrap(angle(h1ModelResponse));
h1TruePhaseResidualRad = h1TruePhase - h1Design * ...
    (h1Design(h1Reliable, :) \ h1TruePhase(h1Reliable));
h1ModelPhaseResidualRad = h1ModelPhase - h1Design * ...
    (h1Design(h1Reliable, :) \ h1ModelPhase(h1Reliable));

% 各路分别归一化到单位 RMS，PSD 只比较 PA 输出恢复的频谱形状。
psdSignals = {y0(idx), yPhase44(idx), yPhaseAmplitude44(idx)};
psdNames = {'真实 PA 输出', '仅相位恢复', '相位+幅度恢复'};
psdDb = zeros(nfft, numel(psdSignals));
for iSignal = 1:numel(psdSignals)
    signal = psdSignals{iSignal};
    [density, psdFrequencyMHz] = pwelch( ...
        signal / rms(signal), options{:});
    psdDb(:, iSignal) = 10 * log10(max(density, realmin));
end

results = struct('h1Order', h1Order, 'h1RippleDb', RpDb, ...
    'h2EdgesMHz', H2edgesMHz, 'h2Orders', h2Orders, ...
    'frequencyMHz', fReliable, 'phaseResidualRad', phaseResidual, ...
    'truePhaseResidualRad', truePhaseResidual, ...
    'magnitudeDb', magnitudeDb, 'trueMagnitudeDb', trueMagnitudeDb, ...
    'h1FrequencyMHz', h1FrequencyMHz, ...
    'h1TrueMagnitudeDb', h1TrueMagnitudeDb, ...
    'h1ModelMagnitudeDb', h1ModelMagnitudeDb, ...
    'h1TruePhaseResidualRad', h1TruePhaseResidualRad, ...
    'h1ModelPhaseResidualRad', h1ModelPhaseResidualRad, ...
    'phaseErrorRmsRad', phaseErrorRmsRad, ...
    'magnitudeErrorRmsDb', magnitudeErrorRmsDb, ...
    'methodNames', {methodNames}, 'nmseDb', nmseDb, ...
    'psdFrequencyMHz', psdFrequencyMHz, 'psdDb', psdDb, ...
    'psdNames', {psdNames}, ...
    'evaluationIndices', [idx(1), idx(end)]);
save(fullfile(resultDir, 'results.mat'), 'results');

%% 可视化：局部校准观测与两种 H2 带宽下的 PA 恢复
fig = plot_calibration_response(results);
exportgraphics(fig, fullfile(resultDir, 'calibration_response.png'), ...
    'Resolution', 180);

fig = figure('Visible', 'off');
bar(nmseDb.');
xticks(1:numel(methodNames));
xticklabels(methodNames);
legend(compose('H2 %d MHz', H2edgesMHz), 'Location', 'best');
ylabel('PA 输出恢复 NMSE (dB)');
title('H1 通带波纹 1 dB：相位与幅度联合校准');
grid on;
exportgraphics(fig, fullfile(resultDir, 'nmse_comparison.png'), ...
    'Resolution', 180);
close(fig);

fig = plot_psd_comparison(results);
exportgraphics(fig, fullfile(resultDir, 'psd_comparison.png'), ...
    'Resolution', 180);

fprintf('H1 %d 阶，波纹 %.1f dB；可靠频段 %.2f 至 %.2f MHz\n', ...
    h1Order, RpDb, min(fReliable), max(fReliable));
fprintf('幅度误差 %.4f dB；残余相位误差 %.5f rad\n', ...
    magnitudeErrorRmsDb, phaseErrorRmsRad);
for iBw = 1:numel(H2edgesMHz)
    fprintf('H2 %d MHz：仅 H2 %.2f；相位 %.2f；相位+幅度 %.2f；参考 %.2f dB\n', ...
        H2edgesMHz(iBw), nmseDb(iBw, :));
end
