%% 输入与反馈相位估计，再经 H2 覆盖恢复 PA 全带输出
% x 为 PA 输入，y0 为实测 PA 输出；H1 模拟未知反馈链，H2 是已知的反馈低通。
% 开关 2 绕过 PA：x_f=H1(x)，由 x,x_f 训练 H1 残余相位。
% 开关 1 接入 PA：实际反馈为 H2*H1(y0)，相位 FIR 仅在模型侧与 H2 组合。
% x,y_f 与 x,x_f 组不使用真实 y0 估计相位；y,y_f 组是离线参考。
clear;
close all;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataDir = fullfile(root, 'datasets', '测试', '2026.6.10', 'test1');
outDir = fullfile(root, 'projects', 'Phase_Calibration', ...
    'results', 'phase_corrected_h2_cover');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

Fs = 1000;                     % MHz；每个采样间隔为 1 ns
BW = 80;                       % H2 正频率通带边缘，MHz
RpDb = 0.1;                    % H1 通带波纹，dB
RsDb = 80;
Dpass2 = 0.0057501127785;
Dstop2 = 1e-4;
K = 7;
M = 5;
S = 0.3;
lamda = 1;
guard = 5000;

d = load(fullfile(dataDir, 'x_v0.6_v1.mat'));
x = d.x(:);
d = load(fullfile(dataDir, 'y0_v0.6_v1.mat'));
y0 = d.y0(:);

%% 观测链路：0.1 dB 椭圆 IIR H1 后接已知 H2
% 实验未逐块模拟 BPF、混频和 ADC，仅以 H1 表示未知反馈频响。
[h1Order, wn] = ellipord(85 / (Fs / 2), 90 / (Fs / 2), RpDb, RsDb);
[b1, a1] = ellip(h1Order, RpDb, RsDb, wn, 'low');
zRaw = filter(b1, a1, y0);

% 对齐只消除整体时延与常相位，保留随频率变化的反馈链相位。
[z, xAligned, ~] = align_circular_signal(zRaw, x);
[h2Order, fo, ao, weights] = firpmord( ...
    [BW, BW + 3] / (Fs / 2), [1 0], [Dpass2, Dstop2]);
h2 = firpm(h2Order, fo, ao, weights, {20});
yFeedback = filter(h2, 1, z);
X = geth_even_1d(xAligned, K, M);

%% 仅在输入激励可靠的频点估计 x 到 z 的残余相位
% S_zx/S_xx 同时含 PA 的线性近似响应与 H1，因而不是纯 H1 测量。
% 只保留输入功率足够、相干度足够的频点，再去掉最佳常相位和线性时延。
idx = (guard + 1):(numel(xAligned) - guard);
segmentLength = 4096;
overlap = segmentLength / 2;
nfft = 8192;
smoothingBins = 5;
window = hann(segmentLength, 'periodic');
options = {window, overlap, nfft, Fs, 'centered'};
[Szx, fMHz] = cpsd(z(idx), xAligned(idx), options{:});
Sxx = pwelch(xAligned(idx), options{:});
Szz = pwelch(z(idx), options{:});
Szx = movmean(Szx, smoothingBins);
Sxx = movmean(Sxx, smoothingBins);
Szz = movmean(Szz, smoothingBins);
coherence = min(1, abs(Szx).^2 ./ (Sxx .* Szz));
reliable = abs(fMHz) <= BW & ...
    10 * log10(Sxx / max(Sxx)) >= -30 & coherence >= 0.8;
assert(nnz(reliable) >= 20, '有效激励内的 x,y_f 相位频点不足。');
fReliable = fMHz(reliable);
omegaReliable = 2 * pi * fReliable / Fs;
phaseObserved = unwrap(angle(Szx(reliable) ./ Sxx(reliable)));
design = [ones(nnz(reliable), 1), omegaReliable];
phaseResidual = phaseObserved - design * (design \ phaseObserved);

%% 开关 2 训练：绕过 PA，测量 x_BL=H1(x)，辨识反馈链残余相位
% x_BL 即图中的校准反馈 x_f；用与 x,y_f 组相同的可靠频点公平比较。
% 训练记录不经过 PA，因此相位估计不混入 PA 的残余线性相位。
xBLRaw = filter(b1, a1, x);
[xBL, ~, ~] = align_circular_signal(xBLRaw, xAligned);
SBLx = movmean(cpsd(xBL(idx), xAligned(idx), options{:}), smoothingBins);
phaseCalibrationObserved = unwrap(angle(SBLx(reliable) ./ Sxx(reliable)));
phaseCalibrationResidual = phaseCalibrationObserved - design * ...
    (design \ phaseCalibrationObserved);

%% 先构造纯相位算子，再与基线相同的 H2 卷积覆盖
% 模型侧采用 +phi_H1，重现 H1 对 PA 输出施加的相移，而不是对反馈求逆。
% 可靠频段边缘的 5 MHz 内将相位修正平滑降为零，未观测频段不作修正。
% 4001 抽头 FIR 的实现时延在随后对齐时消除。
phaseTapCount = 4001;
[hPhase, phaseMagnitudeRippleDb] = create_phase_only_fir( ...
    fReliable, phaseResidual, Fs, nfft, phaseTapCount, 5);
hModel = conv(hPhase, h2(:));
[hPhaseCalibration, phaseCalibrationMagnitudeDeviationDb] = create_phase_only_fir( ...
    fReliable, phaseCalibrationResidual, Fs, nfft, phaseTapCount, 5);
hModelCalibration = conv(hPhaseCalibration, h2(:));

% 检查截断 FIR 是否引入明显幅度偏差。
fprintf('H1 %d 阶；可靠相位频段 %.2f 至 %.2f MHz；相位算子幅度偏差最大 %.3f dB\n', ...
    h1Order, min(fReliable), max(fReliable), phaseMagnitudeRippleDb);

%% LD-PLS 对照：基线与局部相位组使用相同观测、基函数和评价区间
% 只改变模型侧算子；基线为 H2，x,y_f 组为 H2*hPhase。
referenceCover = filter(h2, 1, xAligned);
[feedbackH2, ~, ~] = align_circular_signal(yFeedback, referenceCover);
[yH2, ~] = get_LD_PLS(K, M, X, xAligned, h2, ...
    feedbackH2, S, lamda);

referencePhase = filter(hModel, 1, xAligned);
[feedbackPhase, ~, ~] = align_circular_signal(yFeedback, referencePhase);
[yPhase, ~] = get_LD_PLS(K, M, X, xAligned, hModel, ...
    feedbackPhase, S, lamda);

%% 开关 1 工作：相位算子只进入模型侧，实际 PA 反馈仍由 H2 滤波
% 训练后的 H2*hPhaseCalibration 作用于模型基函数；测量反馈仍为 H2*z。
referenceCalibration = filter(hModelCalibration, 1, xAligned);
[feedbackCalibration, ~, ~] = align_circular_signal( ...
    yFeedback, referenceCalibration);
[yCalibration, ~] = get_LD_PLS(K, M, X, xAligned, ...
    hModelCalibration, feedbackCalibration, S, lamda);

%% H2 窄通带对照：将通带边缘从 80 MHz 缩至 44 MHz
% 相位训练频点和 FIR 保持不变，仅更换实际反馈与模型侧共同使用的 H2。
% 阻带边缘为 47 MHz，已处于局部相位 FIR 的边缘过渡附近。
narrowH2MHz = 44;
[h2EffectiveOrder, foEffective, aoEffective, weightsEffective] = ...
    firpmord([narrowH2MHz, narrowH2MHz + 3] / (Fs / 2), ...
    [1 0], [Dpass2, Dstop2]);
h2Effective = firpm(h2EffectiveOrder, foEffective, aoEffective, ...
    weightsEffective, {20});
yFeedbackEffective = filter(h2Effective, 1, z);

referenceEffective = filter(h2Effective, 1, xAligned);
[feedbackEffective, ~, ~] = align_circular_signal( ...
    yFeedbackEffective, referenceEffective);
[yEffective, ~] = get_LD_PLS(K, M, X, xAligned, h2Effective, ...
    feedbackEffective, S, lamda);

hModelEffective = conv(hPhase, h2Effective(:));
referenceEffectivePhase = filter(hModelEffective, 1, xAligned);
[feedbackEffectivePhase, ~, ~] = align_circular_signal( ...
    yFeedbackEffective, referenceEffectivePhase);
[yEffectivePhase, ~] = get_LD_PLS(K, M, X, xAligned, hModelEffective, ...
    feedbackEffectivePhase, S, lamda);

hModelCalibrationEffective = conv(hPhaseCalibration, h2Effective(:));
referenceCalibrationEffective = filter(hModelCalibrationEffective, 1, xAligned);
[feedbackCalibrationEffective, ~, ~] = align_circular_signal( ...
    yFeedbackEffective, referenceCalibrationEffective);
[yCalibrationEffective, ~] = get_LD_PLS(K, M, X, xAligned, ...
    hModelCalibrationEffective, feedbackCalibrationEffective, S, lamda);

%% 离线参考：拟合完整算子，仅提取其相位，再按相同方式经 H2 覆盖
% 真实 y0 只在这里用于拟合 H1 参考相位，不能作为在线校准可用的信息。
% 局部参考沿用同一可靠频段、边缘过渡宽度和抽头数。
fullOperatorLength = 2048;
fullOperatorEta = 1e-6;
[hFull, fullOperatorFitDb] = estimate_full_fir_operator( ...
    y0, zRaw, fullOperatorLength, fullOperatorEta, idx);
fullPhaseObserved = unwrap(angle(freqz(hFull, 1, fReliable, Fs)));
fullPhaseResidual = fullPhaseObserved - design * (design \ fullPhaseObserved);
phaseEstimateDifferenceRmsRad = rms(phaseResidual - fullPhaseResidual);
phaseCalibrationDifferenceRmsRad = rms(phaseCalibrationResidual - fullPhaseResidual);
[hPhaseFull, fullPhaseMagnitudeDeviationDb] = create_phase_only_fir( ...
    fReliable, fullPhaseResidual, Fs, nfft, phaseTapCount, 5);
hModelFull = conv(hPhaseFull, h2(:));

% 离线相位估计保持有效激励频段不变，只将物理反馈和覆盖换为窄 H2。
hModelFullEffective = conv(hPhaseFull, h2Effective(:));
referenceFullEffective = filter(hModelFullEffective, 1, xAligned);
[feedbackFullEffective, ~, ~] = align_circular_signal( ...
    yFeedbackEffective, referenceFullEffective);
[yFullEffective, ~] = get_LD_PLS(K, M, X, xAligned, hModelFullEffective, ...
    feedbackFullEffective, S, lamda);

% 同一离线算子扩展到 H1 的 ±85 MHz 通带，80–85 MHz 平滑归零。
% 这组带外相位来自 y0,y_f 离线拟合，并非 x 的激励测到了 ±80 MHz。
wideBand = abs(fMHz) <= 85;
fWide = fMHz(wideBand);
wideReliable = reliable(wideBand);
wideDesign = [ones(nnz(wideBand), 1), 2 * pi * fWide / Fs];
widePhaseObserved = unwrap(angle(freqz(hFull, 1, fWide, Fs)));
widePhaseResidual = widePhaseObserved - wideDesign * ...
    (wideDesign(wideReliable, :) \ widePhaseObserved(wideReliable));
[hPhaseWide, widePhaseMagnitudeDeviationDb] = create_phase_only_fir( ...
    fWide, widePhaseResidual, Fs, nfft, phaseTapCount, 5);
hModelWide = conv(hPhaseWide, h2(:));

% 局部相位 FIR 在可靠频段边缘归零，量化其与宽频相位 FIR 的通带边缘差异。
fEdge = (narrowH2MHz - 5:0.1:narrowH2MHz).';
localEdge = freqz(hPhaseFull, 1, fEdge, Fs);
wideEdge = freqz(hPhaseWide, 1, fEdge, Fs);
phaseEdgeDifferenceRmsDeg = rms(angle(localEdge .* conj(wideEdge))) * 180 / pi;

% 宽频离线相位保持不变，仅把物理反馈与模型覆盖换为窄 H2。
hModelWideEffective = conv(hPhaseWide, h2Effective(:));
referenceWideEffective = filter(hModelWideEffective, 1, xAligned);
[feedbackWideEffective, ~, ~] = align_circular_signal( ...
    yFeedbackEffective, referenceWideEffective);
[yWideEffective, ~] = get_LD_PLS(K, M, X, xAligned, hModelWideEffective, ...
    feedbackWideEffective, S, lamda);

fullBand = abs(fMHz) <= BW;
wideMagnitudeH2Db = max(abs(20 * log10(abs( ...
    freqz(hPhaseWide, 1, fMHz(fullBand), Fs)))));
omegaFull = 2 * pi * fMHz(fullBand) / Fs;
fullDesign = [ones(nnz(fullBand), 1), omegaFull];
phaseFull = unwrap(angle(freqz(hFull, 1, fMHz(fullBand), Fs)));
phaseTrueFull = unwrap(angle(freqz(b1, a1, fMHz(fullBand), Fs)));
phaseFull = phaseFull - fullDesign * (fullDesign \ phaseFull);
phaseTrueFull = phaseTrueFull - fullDesign * (fullDesign \ phaseTrueFull);
fullOperatorPhaseErrorRad = rms(phaseFull - phaseTrueFull);
referenceFull = filter(hModelFull, 1, xAligned);
[feedbackFull, ~, ~] = align_circular_signal(yFeedback, referenceFull);
[yFull, ~] = get_LD_PLS(K, M, X, xAligned, hModelFull, ...
    feedbackFull, S, lamda);
referenceWide = filter(hModelWide, 1, xAligned);
[feedbackWide, ~, ~] = align_circular_signal(yFeedback, referenceWide);
[yWide, ~] = get_LD_PLS(K, M, X, xAligned, hModelWide, ...
    feedbackWide, S, lamda);

%% 全带输出评分与结果保存
% x,y_f 与开关 2 组仅在此处使用 y0 评分；y,y_f 组还用它辨识离线算子。
% get_NMSE 分别归一化两路信号 RMS，结果不直接衡量绝对增益误差。
[yH2Aligned, ~, ~] = align_circular_signal(yH2, y0);
[yPhaseAligned, ~, ~] = align_circular_signal(yPhase, y0);
[yCalibrationAligned, ~, ~] = align_circular_signal(yCalibration, y0);
[yFullAligned, ~, ~] = align_circular_signal(yFull, y0);
[yWideAligned, ~, ~] = align_circular_signal(yWide, y0);
[yEffectiveAligned, ~, ~] = align_circular_signal(yEffective, y0);
[yEffectivePhaseAligned, ~, ~] = align_circular_signal(yEffectivePhase, y0);
[yCalibrationEffectiveAligned, ~, ~] = align_circular_signal( ...
    yCalibrationEffective, y0);
[yFullEffectiveAligned, ~, ~] = align_circular_signal(yFullEffective, y0);
[yWideEffectiveAligned, ~, ~] = align_circular_signal(yWideEffective, y0);
nmseH2Db = get_NMSE(yH2Aligned(idx), y0(idx));
nmsePhaseDb = get_NMSE(yPhaseAligned(idx), y0(idx));
nmseCalibrationDb = get_NMSE(yCalibrationAligned(idx), y0(idx));
nmseFullDb = get_NMSE(yFullAligned(idx), y0(idx));
nmseWideDb = get_NMSE(yWideAligned(idx), y0(idx));
nmseEffectiveDb = get_NMSE(yEffectiveAligned(idx), y0(idx));
nmseEffectivePhaseDb = get_NMSE(yEffectivePhaseAligned(idx), y0(idx));
nmseCalibrationEffectiveDb = get_NMSE(yCalibrationEffectiveAligned(idx), y0(idx));
nmseFullEffectiveDb = get_NMSE(yFullEffectiveAligned(idx), y0(idx));
nmseWideEffectiveDb = get_NMSE(yWideEffectiveAligned(idx), y0(idx));

results = struct( ...
     'h1Order', h1Order, 'h1RippleDb', RpDb, ...
    'h2Order', h2Order, ...
    'narrowH2EdgeMHz', narrowH2MHz, ...
    'effectiveH2Order', h2EffectiveOrder, ...
    'phaseFrequencyMHz', fReliable, ...
    'estimatedResidualPhaseRad', phaseResidual, ...
    'switch2CalibrationResidualPhaseRad', phaseCalibrationResidual, ...
    'offlineResidualPhaseRad', fullPhaseResidual, ...
    'offlineWideFrequencyMHz', fWide, ...
    'offlineWideResidualPhaseRad', widePhaseResidual, ...
    'phaseEstimateDifferenceRmsRad', phaseEstimateDifferenceRmsRad, ...
    'switch2PhaseDifferenceRmsRad', phaseCalibrationDifferenceRmsRad, ...
    'phaseOperatorMagnitudeDeviationDb', phaseMagnitudeRippleDb, ...
    'switch2PhaseMagnitudeDeviationDb', phaseCalibrationMagnitudeDeviationDb, ...
    'offlinePhaseMagnitudeDeviationDb', fullPhaseMagnitudeDeviationDb, ...
    'offlineWidePhaseMagnitudeDeviationDb', widePhaseMagnitudeDeviationDb, ...
    'offlineWidePhaseMagnitudeH2Db', wideMagnitudeH2Db, ...
    'phaseEdgeDifferenceRmsDeg', phaseEdgeDifferenceRmsDeg, ...
    'nmseH2Db', nmseH2Db, 'nmsePhaseH2Db', nmsePhaseDb, ...
    'nmseSwitch2PhaseH2Db', nmseCalibrationDb, ...
    'nmseOfflinePhaseH2Db', nmseFullDb, ...
    'nmseOfflineWidePhaseH2Db', nmseWideDb, ...
    'nmseEffectiveH2Db', nmseEffectiveDb, ...
    'nmseEffectivePhaseH2Db', nmseEffectivePhaseDb, ...
    'nmseSwitch2PhaseEffectiveH2Db', nmseCalibrationEffectiveDb, ...
    'nmseOfflineEffectivePhaseH2Db', nmseFullEffectiveDb, ...
    'nmseOfflineWideEffectivePhaseH2Db', nmseWideEffectiveDb, ...
    'fullOperatorFitDb', fullOperatorFitDb, ...
    'fullOperatorPhaseErrorRad', fullOperatorPhaseErrorRad, ...
    'fullOperatorLength', fullOperatorLength, ...
    'nmseChangeDb', nmsePhaseDb - nmseH2Db, ...
    'evaluationIndices', [idx(1), idx(end)]);
save(fullfile(outDir, 'results.mat'), 'results');

%% 可视化：相位来源，以及相同反馈条件下的 NMSE 对照
% 颜色始终对应相位来源；相位来源相同的两根柱只改变 H2 通带边缘。
fig = figure('Name', 'Figure 1：残余相位', 'NumberTitle', 'off');
fig.Position(3) = 1100;
movegui(fig, 'northwest');
plot(fWide, widePhaseResidual * 180 / pi, 'k-', 'LineWidth', 1.2);
hold on;
plot(fReliable, phaseResidual * 180 / pi, 'b-', 'LineWidth', 1.2);
plot(fReliable, phaseCalibrationResidual * 180 / pi, ...
    'Color', [0.46 0.28 0.66], 'LineWidth', 1.4);
plot(fReliable, fullPhaseResidual * 180 / pi, 'r-', 'LineWidth', 1.2);
grid on;
xlabel('频率 (MHz)');
ylabel('残余相位 (度)');
title('有效激励与 H1 通带内的残余相位修正');
legend('y,y_f：H1 通带 ±85 MHz', ...
    'x,y_f：±49 MHz', ...
    'x,x_f：±49 MHz', ...
    'y,y_f：±49 MHz', 'Location', 'best');

colorBaseline = [0.35 0.35 0.35];
colorInputFeedback = [0.08 0.43 0.73];
colorSwitch2 = [0.46 0.28 0.66];
colorOfflineLocal = [0.00 0.55 0.47];
colorOfflineWide = [0.83 0.36 0.13];

fig = figure('Name', 'Figure 2：开关 2 相位训练与九组 NMSE', 'NumberTitle', 'off');
fig.Position(3:4) = [1200, 700];
movegui(fig, 'southwest');
nmseValues = [nmseH2Db, nmsePhaseDb, nmseCalibrationDb, nmseFullDb, nmseWideDb, ...
    nmseEffectivePhaseDb, nmseCalibrationEffectiveDb, ...
    nmseFullEffectiveDb, nmseWideEffectiveDb];
bars = barh(1:9, nmseValues, 0.68, 'FaceColor', 'flat');
bars.CData = [colorBaseline; colorInputFeedback; colorSwitch2; ...
    colorOfflineLocal; colorOfflineWide; colorInputFeedback; ...
    colorSwitch2; colorOfflineLocal; colorOfflineWide];
yticks(1:9);
yticklabels({'仅 H2 80 MHz', 'x,y_f | 相位 ±49 | H2 80', ...
    'x,x_f | 相位 ±49 | H2 80', ...
    'y,y_f | 相位 ±49 | H2 80', 'y,y_f | 相位 ±80 | H2 80', ...
    sprintf('x,y_f | 相位 ±49 | H2 %d', narrowH2MHz), ...
    sprintf('x,x_f | 相位 ±49 | H2 %d', narrowH2MHz), ...
    sprintf('y,y_f | 相位 ±49 | H2 %d', narrowH2MHz), ...
    sprintf('y,y_f | 相位 ±80 | H2 %d', narrowH2MHz)});
xlim([-42, 0]);
xticks(-40:10:0);
xlabel('PA 全带恢复 NMSE (dB，越低越好)');
title('开关 2 训练相位后用于模型侧：九组恢复结果');
grid on;
ax = gca;
ax.YDir = 'reverse';
ax.Position = [0.33, 0.12, 0.63, 0.80];
for k = 1:numel(nmseValues)
    text(nmseValues(k) + 0.5, k, sprintf('%.2f dB', nmseValues(k)), ...
        'Color', 'white', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
keyBaseline = patch(NaN, NaN, colorBaseline, 'EdgeColor', 'none');
keyInput = patch(NaN, NaN, colorInputFeedback, 'EdgeColor', 'none');
keySwitch2 = patch(NaN, NaN, colorSwitch2, 'EdgeColor', 'none');
keyLocal = patch(NaN, NaN, colorOfflineLocal, 'EdgeColor', 'none');
keyWide = patch(NaN, NaN, colorOfflineWide, 'EdgeColor', 'none');
legend([keyBaseline, keyInput, keySwitch2, keyLocal, keyWide], ...
    {'仅 H2', 'x,y_f 估计', 'x,x_f', 'y,y_f 局部参考', ...
    'y,y_f 宽频参考'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');

fig = figure('Name', 'Figure 3：开关 2 与 H2 带宽对比', 'NumberTitle', 'off');
fig.Position(3:4) = [1200, 650];
movegui(fig, 'southeast');
nmseValues = [nmseEffectivePhaseDb, nmsePhaseDb, ...
    nmseCalibrationEffectiveDb, nmseCalibrationDb, ...
    nmseFullEffectiveDb, nmseFullDb, ...
    nmseWideEffectiveDb, nmseWideDb];
bars = barh(1:8, nmseValues, 0.68, 'FaceColor', 'flat');
bars.CData = [colorInputFeedback; colorInputFeedback; ...
    colorSwitch2; colorSwitch2; ...
    colorOfflineLocal; colorOfflineLocal; ...
    colorOfflineWide; colorOfflineWide];
yticks(1:8);
yticklabels({sprintf('x,y_f | 相位 ±49 | H2 %d', narrowH2MHz), ...
    'x,y_f | 相位 ±49 | H2 80', ...
    sprintf('x,x_f | 相位 ±49 | H2 %d', narrowH2MHz), ...
    'x,x_f | 相位 ±49 | H2 80', ...
    sprintf('y,y_f | 相位 ±49 | H2 %d', narrowH2MHz), ...
    'y,y_f | 相位 ±49 | H2 80', ...
    sprintf('y,y_f | 相位 ±80 | H2 %d', narrowH2MHz), ...
    'y,y_f | 相位 ±80 | H2 80'});
xlim([-42, 0]);
xticks(-40:10:0);
xlabel('PA 全带恢复 NMSE (dB，越低越好)');
title('开关 2 相位训练与 H2 带宽对比');
grid on;
ax = gca;
ax.YDir = 'reverse';
ax.Position = [0.33, 0.12, 0.63, 0.80];
for k = 1:numel(nmseValues)
    text(nmseValues(k) + 0.5, k, sprintf('%.2f dB', nmseValues(k)), ...
        'Color', 'white', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
keyInput = patch(NaN, NaN, colorInputFeedback, 'EdgeColor', 'none');
keySwitch2 = patch(NaN, NaN, colorSwitch2, 'EdgeColor', 'none');
keyLocal = patch(NaN, NaN, colorOfflineLocal, 'EdgeColor', 'none');
keyWide = patch(NaN, NaN, colorOfflineWide, 'EdgeColor', 'none');
legend([keyInput, keySwitch2, keyLocal, keyWide], ...
    {'x,y_f 估计', 'x,x_f', 'y,y_f 局部参考', ...
    'y,y_f 宽频参考'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
figure(fig);
drawnow;

fprintf('仅 H2 覆盖 NMSE：%.2f dB\n', nmseH2Db);
fprintf('x,y_f 有效激励内相位修正 + H2：%.2f dB；变化 %.2f dB\n', ...
    nmsePhaseDb, results.nmseChangeDb);
fprintf('y,y_f 有效激励内相位修正 + H2：%.2f dB；算子拟合 %.2f dB\n', ...
    nmseFullDb, fullOperatorFitDb);
fprintf('y,y_f H1 通带相位修正 + H2：%.2f dB；比 ±49 MHz 改善 %.2f dB\n', ...
    nmseWideDb, nmseFullDb - nmseWideDb);
fprintf('H2 %d MHz：x,y_f 局部相位修正 %.2f dB；仅 H2 %.2f dB\n', ...
    narrowH2MHz, nmseEffectivePhaseDb, nmseEffectiveDb);
fprintf('开关 2 训练相位，模型侧补偿：H2 %d MHz %.2f dB；H2 %d MHz %.2f dB\n', ...
    narrowH2MHz, nmseCalibrationEffectiveDb, BW, nmseCalibrationDb);
fprintf('H2 %d MHz：y,y_f 局部相位修正 %.2f dB\n', ...
    narrowH2MHz, nmseFullEffectiveDb);
fprintf('H2 %d MHz：y,y_f 宽频相位修正 %.2f dB\n', ...
    narrowH2MHz, nmseWideEffectiveDb);
fprintf('%d–%d MHz 局部与宽频相位 FIR 的差异 RMS：%.2f 度\n', ...
    narrowH2MHz - 5, narrowH2MHz, phaseEdgeDifferenceRmsDeg);
fprintf('两种相位估计相差 RMS %.5f rad；恢复 NMSE 相差 %.2f dB\n', ...
    phaseEstimateDifferenceRmsRad, nmsePhaseDb - nmseFullDb);
fprintf('离线完整算子在 H2 通带内的残余相位误差：%.5f rad\n', ...
    fullOperatorPhaseErrorRad);
fprintf('x,y_f 相位仅在有效激励频段估计；±85 MHz 组使用 y,y_f 的参考算子。\n');

%% 将本次运行的参数和数值写入可追溯报告
reportPath = fullfile(outDir, 'report.md');
fid = fopen(reportPath, 'w');
fprintf(fid, '# 有效激励与 H1 通带相位修正的 H2 覆盖对照\n\n');
fprintf(fid, '实测 PA 输入输出 38400 点，Fs=1000 MHz。H1 为 %d 阶椭圆 IIR，通带/阻带边缘 85/90 MHz，通带波纹 %.1f dB。H2 通带边缘 %d MHz，阶数 %d。\n\n', ...
    h1Order, RpDb, BW, h2Order);
fprintf(fid, '令 y_f=H1(y)。开关 2 绕过 PA，得到 x_BL=H1(x)，先由 x,x_BL 训练反馈链残余相位；切回开关 1 后，实际 PA 反馈保持 H2*H1(y)，训练得到的相位 FIR 仅在模型侧与 H2 卷积覆盖。另由 x,y_f 互谱估计相位。离线参考用真实 y 与原始 y_f 拟合 %d 抽头完整复数 FIR，但只提取相位。局部三组修正约 ±49 MHz；宽频组使用同一 y,y_f 算子修正至 ±80 MHz，并在 80–85 MHz 过渡归零。各相位组使用相同的 %d 抽头纯相位 FIR 构造方法。前五组使用同一份 H2 80 MHz 后反馈；后四组使用同一份 H2 %d MHz 后反馈。九组使用相同 LD-PLS 参数，并在样本 %d–%d 上评分。输入 x 的可靠激励始终约为 ±49 MHz；宽频组扩展的是离线相位修正范围，并非输入激励范围。开关 2 是 PA 旁路校准，不是实际功率回退数据。\n\n', ...
    fullOperatorLength, phaseTapCount, narrowH2MHz, ...
    idx(1), idx(end));
fprintf(fid, '| 方法 | PA 全带恢复 NMSE |\n|---|---:|\n');
fprintf(fid, '| 仅 H2 80 MHz 覆盖 | %.2f dB |\n', nmseH2Db);
fprintf(fid, '| x,y_f：相位约 ±49 MHz + H2 80 MHz 覆盖 | %.2f dB |\n', nmsePhaseDb);
fprintf(fid, '| 开关 2：x,x_BL 训练相位 + H2 80 MHz 模型侧覆盖 | %.2f dB |\n', nmseCalibrationDb);
fprintf(fid, '| y,y_f：相位约 ±49 MHz + H2 80 MHz 覆盖 | %.2f dB |\n', nmseFullDb);
fprintf(fid, '| y,y_f：相位修正至 ±80 MHz + H2 80 MHz 覆盖 | %.2f dB |\n', nmseWideDb);
fprintf(fid, '| x,y_f：相位约 ±49 MHz + H2 %d MHz 覆盖 | %.2f dB |\n', narrowH2MHz, nmseEffectivePhaseDb);
fprintf(fid, '| 开关 2：x,x_BL 训练相位 + H2 %d MHz 模型侧覆盖 | %.2f dB |\n', narrowH2MHz, nmseCalibrationEffectiveDb);
fprintf(fid, '| y,y_f：相位约 ±49 MHz + H2 %d MHz 覆盖 | %.2f dB |\n', narrowH2MHz, nmseFullEffectiveDb);
fprintf(fid, '| y,y_f：相位修正至 ±80 MHz + H2 %d MHz 覆盖 | %.2f dB |\n\n', narrowH2MHz, nmseWideEffectiveDb);
fprintf(fid, '## 相同输入记录下的相位修正范围与 H2 带宽对比\n\n');
fprintf(fid, '同一条 x、y 实测记录保持不变，x 的可靠激励边缘实测为 %.2f MHz；窄 H2 通带/阻带边缘为 %d/%d MHz，阶数 %d。物理反馈和模型侧覆盖均使用对应的 H2；宽频 y,y_f 组把离线相位修正扩展至 ±80 MHz，输入激励没有变宽。所有组保持 LD-PLS 参数和评分样本一致，不增加过拟合处理。\n\n', ...
    max(abs(fReliable)), narrowH2MHz, narrowH2MHz + 3, h2EffectiveOrder);
fprintf(fid, '| 相位来源与修正范围 | H2 通带 | PA 全带恢复 NMSE |\n');
fprintf(fid, '|---|---:|---:|\n');
fprintf(fid, '| x,y_f；相位约 ±49 MHz | ±%d MHz | %.2f dB |\n', ...
    narrowH2MHz, nmseEffectivePhaseDb);
fprintf(fid, '| x,y_f；相位约 ±49 MHz | ±%d MHz | %.2f dB |\n', ...
    BW, nmsePhaseDb);
fprintf(fid, '| 开关 2；x,x_BL 训练相位约 ±49 MHz | ±%d MHz | %.2f dB |\n', ...
    narrowH2MHz, nmseCalibrationEffectiveDb);
fprintf(fid, '| 开关 2；x,x_BL 训练相位约 ±49 MHz | ±%d MHz | %.2f dB |\n', ...
    BW, nmseCalibrationDb);
fprintf(fid, '| y,y_f；相位约 ±49 MHz | ±%d MHz | %.2f dB |\n', ...
    narrowH2MHz, nmseFullEffectiveDb);
fprintf(fid, '| y,y_f；相位约 ±49 MHz | ±%d MHz | %.2f dB |\n', ...
    BW, nmseFullDb);
fprintf(fid, '| y,y_f；相位至 ±80 MHz | ±%d MHz | %.2f dB |\n', ...
    narrowH2MHz, nmseWideEffectiveDb);
fprintf(fid, '| y,y_f；相位至 ±80 MHz | ±%d MHz | %.2f dB |\n\n', ...
    BW, nmseWideDb);
fprintf(fid, '与 H2 80 MHz 相比，H2 %d MHz 下相同局部相位修正的 NMSE 变化 %.2f dB；负值表示改善。仅 H2 覆盖的诊断值分别为 %.2f dB（%d MHz）和 %.2f dB（80 MHz）。\n\n', ...
    narrowH2MHz, nmseEffectivePhaseDb - nmsePhaseDb, ...
    nmseEffectiveDb, narrowH2MHz, nmseH2Db);
fprintf(fid, 'H2 %d MHz 下，y,y_f 局部相位修正相对 x,y_f 局部相位修正的 NMSE 变化 %.2f dB。\n\n', ...
    narrowH2MHz, nmseFullEffectiveDb - nmseEffectivePhaseDb);
fprintf(fid, '固定 H2 %d MHz 时，y,y_f 宽频相位修正相对约 ±49 MHz 相位修正的 NMSE 变化 %.2f dB；固定宽频相位修正时，H2 %d MHz 相对 H2 80 MHz 的 NMSE 变化 %.2f dB。负值表示改善。\n\n', ...
    narrowH2MHz, nmseWideEffectiveDb - nmseFullEffectiveDb, ...
    narrowH2MHz, nmseWideEffectiveDb - nmseWideDb);
fprintf(fid, '局部相位 FIR 在可靠激励边缘的 5 MHz 内逐渐归零；宽频相位 FIR 则继续修正 H2 %d MHz 通带边缘。两者在 %d–%d MHz 的相位差 RMS 为 %.2f 度。\n\n', ...
    narrowH2MHz, narrowH2MHz - 5, narrowH2MHz, ...
    phaseEdgeDifferenceRmsDeg);
fprintf(fid, '固定 H2 80 MHz 时，y,y_f 相位修正由约 ±49 MHz 扩展至 ±80 MHz 后，NMSE 改善 %.2f dB。x,y_f 局部相位修正比 y,y_f 局部相位修正差 %.2f dB。\n\n', ...
    nmseFullDb - nmseWideDb, nmsePhaseDb - nmseFullDb);
fprintf(fid, 'x,y_f 相位修正使 NMSE 改善 %.2f dB；三个相位 FIR 在各自设计频段的最大幅度偏差分别为 %.3f、%.3f 和 %.3f dB，其中宽频段相位 FIR 在 H2 的 ±80 MHz 内为 %.3f dB。\n\n', ...
    nmseH2Db - nmsePhaseDb, phaseMagnitudeRippleDb, ...
    fullPhaseMagnitudeDeviationDb, widePhaseMagnitudeDeviationDb, ...
    wideMagnitudeH2Db);
fprintf(fid, 'x,y_f 与 y,y_f 在相同有效激励频点的残余相位差 RMS 为 %.5f rad；y,y_f 相位修正比 x,y_f 相位修正的恢复 NMSE 低 %.2f dB。\n\n', ...
    phaseEstimateDifferenceRmsRad, nmsePhaseDb - nmseFullDb);
fprintf(fid, '开关 2 训练相位与 y,y_f 离线相位在相同频点的残余相位差 RMS 为 %.6f rad。H2 %d MHz 下，开关 2 模型侧相位补偿比 x,y_f 相位估计的 NMSE 改善 %.2f dB，与 y,y_f 局部离线参考相差 %.2f dB。\n\n', ...
    phaseCalibrationDifferenceRmsRad, narrowH2MHz, ...
    nmseEffectivePhaseDb - nmseCalibrationEffectiveDb, ...
    nmseCalibrationEffectiveDb - nmseFullEffectiveDb);
fprintf(fid, 'y,y_f 完整算子对 y_f 的时域拟合 NMSE 为 %.2f dB，H2 通带内相对真实 H1 的残余相位误差 RMS 为 %.5f rad；拟合长度 %d 抽头、岭参数 %.0e。y,y_f 各组只用该算子的相位，幅度不进入模型侧覆盖；算子拟合使用真实 y，不能解释为仅凭 x,y_f 可实现的盲校准。\n\n', ...
    fullOperatorFitDb, fullOperatorPhaseErrorRad, ...
    fullOperatorLength, fullOperatorEta);
fprintf(fid, 'x,y_f 相位仅在 %.2f 至 %.2f MHz 的有效激励频段估计，频段边缘外的修正逐渐归零；±85 MHz 组的相位来自真实 y,y_f 拟合算子，并非 x,y_f 的全带估计。\n', ...
    min(fReliable), max(fReliable));
fclose(fid);
