function results = verify_input_feedback_phase()
%VERIFY_INPUT_FEEDBACK_PHASE 验证 x,y_f 的 BLA 是否能近似 H1 残余相位。
%   估计器只用 x,y_f；真实 y 和 H1 仅用于离线误差归因。

    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    dataDir = fullfile(root, 'datasets', '测试', '2026.6.10', 'test1');
    outDir = fullfile(root, 'projects', 'Phase_Calibration', ...
        'results', 'input_feedback_phase');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    d = load(fullfile(dataDir, 'x_v0.6_v1.mat'));
    x = d.x(:);
    d = load(fullfile(dataDir, 'y0_v0.6_v1.mat'));
    y = d.y0(:);

    Fs = 1000;                 % MHz；每个采样间隔为 1 ns
    passEdgeMHz = 85;
    stopEdgeMHz = 90;
    RpDb = 0.1;
    RsDb = 80;
    phaseGateRad = 0.10;
    guard = 5000;              % 排除滤波启动与循环对齐的边缘
    segmentLength = 4096;
    overlap = segmentLength / 2;
    nfft = 8192;
    smoothingBins = 5;
    inputPowerFloorDb = -30;
    coherenceFloor = 0.8;

    [order, wn] = ellipord(passEdgeMHz / (Fs / 2), ...
        stopEdgeMHz / (Fs / 2), RpDb, RsDb);
    [b1, a1] = ellip(order, RpDb, RsDb, wn, 'low');
    z = filter(b1, a1, y);

    % 时间对齐与频谱估计只使用 x 和反馈 z（即 y_f）。
    [zAligned, ~] = align_circular_signal(z, x);
    idx = (guard + 1):(numel(x) - guard);
    xValid = x(idx);
    zValid = zAligned(idx);
    window = hann(segmentLength, 'periodic');
    options = {window, overlap, nfft, Fs, 'centered'};
    [Szx, fMHz] = cpsd(zValid, xValid, options{:});
    Sxx = pwelch(xValid, options{:});
    Szz = pwelch(zValid, options{:});
    Szx = movmean(Szx, smoothingBins);
    Sxx = movmean(Sxx, smoothingBins);
    Szz = movmean(Szz, smoothingBins);
    Gzx = Szx ./ Sxx;
    coherence = min(1, abs(Szx).^2 ./ (Sxx .* Szz));

    inputReliable = 10 * log10(Sxx / max(Sxx)) >= inputPowerFloorDb;
    feedbackReliable = coherence >= coherenceFloor;
    estimatorMask = inputReliable & feedbackReliable;

    % 完成 x,y_f 估计后，才引入真实 H1 进行离线评价。
    H1 = freqz(b1, a1, fMHz, Fs);
    passband = abs(fMHz) <= passEdgeMHz & ...
        abs(H1) >= 10^(-1 / 20);
    mask = estimatorMask & passband;
    if nnz(mask) < 20
        error('verify_input_feedback_phase:InsufficientBins', ...
            '预定输入功率及相干度门限下可靠频点不足。');
    end

    % 互谱方向已经用已知相移信号验证，保持与主程序一致。
    omega = 2 * pi * fMHz(mask) / Fs;
    phiEstimated = residualPhase(Gzx(mask), omega);
    phiTrue = residualPhase(H1(mask), omega);
    errorPhase = phiEstimated - phiTrue;
    errorRmsRad = sqrt(mean(errorPhase.^2));

    % 离线诊断：真实 PA 输出不进入相位估计器。
    [yAligned, ~] = align_circular_signal(y, x);
    yValid = yAligned(idx);
    Syx = movmean(cpsd(yValid, xValid, options{:}), smoothingBins);
    Szy = movmean(cpsd(zValid, yValid, options{:}), smoothingBins);
    Syy = movmean(pwelch(yValid, options{:}), smoothingBins);
    Gyx = Syx ./ Sxx;
    Gzy = Szy ./ Syy;
    phiPa = residualPhase(Gyx(mask), omega);
    phiYz = residualPhase(Gzy(mask), omega);
    spectralError = errorPhase - phiPa;
    spectralErrorYz = phiYz - phiTrue;
    paRmsRad = sqrt(mean(phiPa.^2));
    spectralRmsRad = sqrt(mean(spectralError.^2));
    yzRmsRad = sqrt(mean(spectralErrorYz.^2));
    crossTerm = 2 * mean(phiPa .* spectralError);

    results = struct( ...
        'h1Order', order, 'h1PassbandRippleDb', RpDb, ...
        'frequencyMHz', fMHz(mask), ...
        'evaluationBinCount', nnz(mask), ...
        'passbandBinCount', nnz(passband), ...
        'passbandCoverage', nnz(mask) / nnz(passband), ...
        'fullPassbandValidated', all(estimatorMask(passband)), ...
        'phasePaRad', phiPa, 'phaseH1TrueRad', phiTrue, ...
        'phaseH1EstimatedRad', phiEstimated, ...
        'phaseErrorRad', errorPhase, ...
        'h1TrueRmsRad', sqrt(mean(phiTrue.^2)), ...
        'phaseErrorRmsRad', errorRmsRad, ...
        'paRmsRad', paRmsRad, ...
        'spectralErrorRmsRad', spectralRmsRad, ...
        'yzDiagnosticRmsRad', yzRmsRad, ...
        'decompositionCrossTermRad2', crossTerm, ...
        'phaseGateRad', phaseGateRad, ...
        'passesObservedBandGate', errorRmsRad <= phaseGateRad);

    makeFigures(results);
    save(fullfile(outDir, 'results.mat'), 'results');
    writeReport(outDir, results, segmentLength, overlap, nfft, ...
        smoothingBins, guard, inputPowerFloorDb, coherenceFloor);

    fprintf('H1 %d 阶，0.1 dB；可靠频率 %.2f 到 %.2f MHz，通带覆盖 %.1f%%\n', ...
        order, min(results.frequencyMHz), max(results.frequencyMHz), ...
        100 * results.passbandCoverage);
    fprintf('PA 残余相位 RMS %.4f rad；H1 估计误差 RMS %.4f rad (%.2f deg)\n', ...
        paRmsRad, errorRmsRad, errorRmsRad * 180 / pi);
    fprintf('频谱估计偏差 RMS %.4f rad；y,z 离线检查 %.4f rad\n', ...
        spectralRmsRad, yzRmsRad);
    fprintf('已观测频段相位门限 %.2f rad：%s；全通带验证：%s\n', ...
        phaseGateRad, string(results.passesObservedBandGate), ...
        string(results.fullPassbandValidated));
end

function phase = residualPhase(response, omega)
% 在同一组频点去除最佳常相位与线性时延。
    unwrapped = unwrap(angle(response));
    design = [ones(numel(omega), 1), omega(:)];
    phase = unwrapped - design * (design \ unwrapped);
end

function makeFigures(results)
    f = results.frequencyMHz;
    figureHandle = figure;
    figureHandle.Position(3:4) = [900, 500];
    movegui(figureHandle, 'center');
    plot(f, results.phasePaRad * 180 / pi, 'b-', 'LineWidth', 1.2);
    grid on;
    xlabel('频率 (MHz)');
    ylabel('残余相位 (度)');
    title('实测 PA 的 BLA 残余相位');
    figureHandle = figure;
    figureHandle.Position(3:4) = [900, 500];
    movegui(figureHandle, 'center');
    plot(f, results.phaseH1TrueRad * 180 / pi, 'k-', 'LineWidth', 1.3);
    hold on;
    plot(f, results.phaseH1EstimatedRad * 180 / pi, 'r-', 'LineWidth', 1.1);
    grid on;
    xlabel('频率 (MHz)');
    ylabel('残余相位 (度)');
    title('H1 真实与 x,y_f 估计的残余相位');
    legend('H1 真值', 'x,y_f 估计', 'Location', 'best');
end

function writeReport(outDir, r, segmentLength, overlap, nfft, ...
        smoothingBins, guard, inputPowerFloorDb, coherenceFloor)
    reportPath = fullfile(outDir, 'report.md');
    fid = fopen(reportPath, 'w');
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '# 输入与带限反馈的 H1 残余相位验证\n\n');
    fprintf(fid, '## 实验设置\n\n');
    fprintf(fid, '实测 PA 输入输出 38400 点，Fs=1000 MHz；H1 为 %d 阶椭圆 IIR，通带/阻带边缘 85/90 MHz，通带波纹 %.1f dB，阻带衰减 80 dB。z=H1(y)，无 H2。\n\n', ...
        r.h1Order, r.h1PassbandRippleDb);
    fprintf(fid, '对齐后首尾各去除 %d 点。Welch 分段 %d 点 Hann 窗、重叠 %d 点、NFFT=%d；互谱和功率谱再作 %d 频点移动平均。\n\n', ...
        guard, segmentLength, overlap, nfft, smoothingBins);
    fprintf(fid, '估计频点预设门限：输入功率相对峰值 ≥ %d dB、x,z 相干度 ≥ %.1f。离线评价再限制于 H1 的 ±85 MHz、幅度 ≥ −1 dB 通带。所有相位在同一频点展开，再分别去除最佳常相位与线性时延。\n\n', ...
        inputPowerFloorDb, coherenceFloor);
    fprintf(fid, '## 实测结果\n\n');
    fprintf(fid, '- 可靠频率：%.2f 至 %.2f MHz；%d/%d 个通带频点，覆盖 %.1f%%。\n', ...
        min(r.frequencyMHz), max(r.frequencyMHz), r.evaluationBinCount, ...
        r.passbandBinCount, 100 * r.passbandCoverage);
    fprintf(fid, '- PA BLA 残余相位 RMS：%.5f rad (%.2f°)。\n', ...
        r.paRmsRad, r.paRmsRad * 180 / pi);
    fprintf(fid, '- H1 真实残余相位 RMS：%.5f rad (%.2f°)。\n', ...
        r.h1TrueRmsRad, r.h1TrueRmsRad * 180 / pi);
    fprintf(fid, '- 已观测频段 H1 估计误差 RMS：%.5f rad (%.2f°)，为该频段真值残余相位 RMS 的 %.1f%%；门限 %.2f rad，%s。\n', ...
        r.phaseErrorRmsRad, r.phaseErrorRmsRad * 180 / pi, ...
        100 * r.phaseErrorRmsRad / r.h1TrueRmsRad, ...
        r.phaseGateRad, passText(r.passesObservedBandGate));
    fprintf(fid, '- 频谱估计偏差 RMS：%.5f rad；y,z 离线滤波器相位检查 RMS：%.5f rad。\n', ...
        r.spectralErrorRmsRad, r.yzDiagnosticRmsRad);
    fprintf(fid, '- PA 与频谱偏差的二倍交叉项：%.6f rad²；两项可相互增强或抵消。\n\n', ...
        r.decompositionCrossTermRad2);
    fprintf(fid, '误差分解采用同一去线性算子：估计误差 = PA BLA 残余相位 + 频谱估计偏差。后者由 x,z 与 x,y 的离线差值得到；y,z 互谱提供独立检查。\n\n');
    if r.paRmsRad > r.spectralErrorRmsRad
        source = 'PA 自身残余相位';
    else
        source = '互谱估计偏差';
    end
    fprintf(fid, '## 判断\n\n');
    fprintf(fid, 'PA 残余相位 RMS 约为已观测频段 H1 真实残余相位的 %.1f%%；按分量 RMS 比较，本次主要误差来源为%s。', ...
        100 * r.paRmsRad / r.h1TrueRmsRad, source);
    if r.passesObservedBandGate
        fprintf(fid, 'x,z 的估计在已观测频段达到 0.10 rad 门限。');
    else
        fprintf(fid, 'x,z 的估计在已观测频段未达到 0.10 rad 门限。');
    end
    fprintf(fid, '现有 x 在通带外侧缺乏可靠激励，完整 ±85 MHz 通带尚未验证；不能据此判定全通带估计成功，也不开展联合校准。需要覆盖整个通带的实测 x,y 记录。\n\n');
    fprintf(fid, '本次结果只覆盖一条实测 PA 记录和一个 H1 设计，可靠频段之外不作相位推断。\n');
end

function label = passText(passed)
    if passed
        label = '达到';
    else
        label = '未达到';
    end
end
