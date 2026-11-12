function results = verify_linear_backoff_phase()
%VERIFY_LINEAR_BACKOFF_PHASE 用 x_f=H1(x) 模拟理想线性区的相位校准。
%   校准相位仅由 x,x_f 得到，再用于真实非线性 PA 的反馈恢复。
%   同时检查线性链路自身的恢复；真实 y 和 H1 相位仅作离线对照。

    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    dataDir = fullfile(root, 'datasets', '测试', '2026.6.10', 'test1');
    outDir = fullfile(root, 'projects', 'Phase_Calibration', ...
        'results', 'linear_backoff_phase');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    d = load(fullfile(dataDir, 'x_v0.6_v1.mat'));
    x = d.x(:);
    d = load(fullfile(dataDir, 'y0_v0.6_v1.mat'));
    y = d.y0(:);
    Fs = 1000;
    K = 7;
    M = 5;
    S = 0.3;
    lamda = 1;
    h2EdgeMHz = 45;
    idx = 5001:(numel(x) - 5000);
    nfft = 8192;
    phaseTapCount = 4001;
    options = {hann(4096, 'periodic'), 2048, nfft, Fs, 'centered'};

    %% 生成两条反馈：线性区 x_f，以及真实 PA 的 y_f
    [h1Order, wn] = ellipord(85 / (Fs / 2), 90 / (Fs / 2), 0.1, 80);
    [b1, a1] = ellip(h1Order, 0.1, 80, wn, 'low');
    xfRaw = filter(b1, a1, x);
    yfRaw = filter(b1, a1, y);
    [xf, xAligned, ~] = align_circular_signal(xfRaw, x);
    [yf, ~, ~] = align_circular_signal(yfRaw, xAligned);
    [h2Order, fo, ao, weights] = firpmord( ...
        [h2EdgeMHz, h2EdgeMHz + 3] / (Fs / 2), ...
        [1 0], [0.0057501127785, 1e-4]);
    h2 = firpm(h2Order, fo, ao, weights, {20});

    %% 仅从 x,x_f 选择可靠频点并估计相位，不传入 y 或 H1
    [Sfx, f] = cpsd(xf(idx), xAligned(idx), options{:});
    Sfx = movmean(Sfx, 5);
    Sxx = movmean(pwelch(xAligned(idx), options{:}), 5);
    Sff = movmean(pwelch(xf(idx), options{:}), 5);
    coherence = min(1, abs(Sfx).^2 ./ (Sxx .* Sff));
    reliable = abs(f) <= 80 & ...
        10 * log10(Sxx / max(Sxx)) >= -30 & coherence >= 0.8;
    fPhase = f(reliable);
    design = [ones(numel(fPhase), 1), 2 * pi * fPhase / Fs];
    phaseLinear = residualPhase(Sfx(reliable) ./ Sxx(reliable), design);

    %% 在同一频点生成非线性 PA、离线算子和 H1 真值对照
    Syfx = movmean(cpsd(yf(idx), xAligned(idx), options{:}), 5);
    phaseNonlinear = residualPhase(Syfx(reliable) ./ Sxx(reliable), design);
    [hOffline, offlineFitDb] = estimate_full_fir_operator( ...
        y, yfRaw, 2048, 1e-6, idx);
    phaseOffline = residualPhase(freqz(hOffline, 1, fPhase, Fs), design);
    phaseTrue = residualPhase(freqz(b1, a1, fPhase, Fs), design);
    phaseList = [phaseNonlinear, phaseLinear, phaseOffline, phaseTrue];
    models = cell(1, 5);
    models{1} = h2;
    magnitudeDeviationDb = zeros(1, 4);
    for k = 1:4
        [hPhase, magnitudeDeviationDb(k)] = create_phase_only_fir( ...
            fPhase, phaseList(:, k), Fs, nfft, phaseTapCount, 5);
        models{k + 1} = conv(hPhase, h2(:));
    end

    %% 同一真实 PA 反馈、H2、LD-PLS 参数及评分区间
    X = geth_even_1d(xAligned, K, M);
    measuredFeedback = filter(h2, 1, yf);
    nmseDb = zeros(1, numel(models));
    for k = 1:numel(models)
        nmseDb(k) = recoverNmse(xAligned, y, measuredFeedback, ...
            models{k}, X, K, M, S, lamda, idx);
    end

    % 线性链路一致性检查：此时被恢复对象就是 x，不含 PA 非线性。
    linearFeedback = filter(h2, 1, xf);
    linearModelIndices = [1, 3, 5];
    linearChainNmseDb = zeros(1, 3);
    for k = 1:3
        linearChainNmseDb(k) = recoverNmse(xAligned, xAligned, linearFeedback, ...
            models{linearModelIndices(k)}, X, K, M, S, lamda, idx);
    end

    labels = {'仅 H2 45 MHz', 'x,y_f | 相位 ±49 | H2 45', ...
        'x,x_f | 相位 ±49 | H2 45', 'y,y_f | 相位 ±49 | H2 45', ...
        '真实 H1 相位 ±49 | H2 45'};
    results = struct('h1Order', h1Order, 'h2EdgeMHz', h2EdgeMHz, ...
        'frequencyMHz', fPhase, 'phaseLinearRad', phaseLinear, ...
        'phaseNonlinearRad', phaseNonlinear, 'phaseOfflineRad', phaseOffline, ...
        'phaseTrueRad', phaseTrue, 'phaseErrorRmsRad', ...
        sqrt(mean((phaseList - phaseTrue).^2, 1)), ...
        'phaseOperatorMagnitudeDeviationDb', magnitudeDeviationDb, ...
        'methodLabels', {labels}, 'nmseDb', nmseDb, ...
        'linearChainNmseDb', linearChainNmseDb, 'offlineOperatorFitDb', offlineFitDb, ...
        'evaluationIndices', [idx(1), idx(end)]);
    save(fullfile(outDir, 'results.mat'), 'results');
    writeReport(fullfile(outDir, 'report.md'), results);
    fprintf('可靠频率 %.2f 至 %.2f MHz；H2 45 MHz\n', min(fPhase), max(fPhase));
    fprintf('x,x_f 相位误差 %.6f rad（%.4f 度）；x,y_f %.6f rad\n', ...
        results.phaseErrorRmsRad(2), results.phaseErrorRmsRad(2) * 180 / pi, ...
        results.phaseErrorRmsRad(1));
    for k = 1:numel(labels)
        fprintf('%s：PA 全带恢复 %.2f dB\n', labels{k}, nmseDb(k));
    end
    fprintf('线性链路恢复：仅 H2 %.2f；x,x_f 相位 %.2f；真实相位 %.2f dB\n', ...
        linearChainNmseDb);

    fig = figure;
    fig.Position(3:4) = [900, 500];
    movegui(fig, 'center');
    plot(fPhase, [phaseTrue, phaseNonlinear, phaseLinear] * 180 / pi, ...
        'LineWidth', 1.3);
    grid on;
    xlabel('频率 (MHz)');
    ylabel('残余相位 (度)');
    title('理想线性区 x,x_f 与非线性 PA x,y_f 的相位估计');
    legend('H1 真值', 'x,y_f', 'x,x_f', 'Location', 'best');

    fig = figure;
    fig.Position(3:4) = [1050, 540];
    movegui(fig, 'center');
    bars = barh(1:numel(labels), nmseDb, 0.65, 'FaceColor', 'flat');
    bars.CData = [0.35 0.35 0.35; 0.08 0.43 0.73; 0.46 0.28 0.66; ...
        0.00 0.55 0.47; 0.83 0.36 0.13];
    yticks(1:numel(labels));
    yticklabels(labels);
    ax = gca;
    ax.YDir = 'reverse';
    ax.Position = [0.35, 0.14, 0.61, 0.77];
    xlim([-42, 0]);
    grid on;
    xlabel('PA 全带恢复 NMSE (dB，越低越好)');
    title('线性区校准相位用于真实 PA 恢复：H2 45 MHz');
    for k = 1:numel(labels)
        text(nmseDb(k) + 0.5, k, sprintf('%.2f dB', nmseDb(k)), ...
            'Color', 'white', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'middle');
    end
end

function phase = residualPhase(response, design)
    phase = unwrap(angle(response));
    phase = phase - design * (design \ phase);
end

function nmseDb = recoverNmse(x, target, feedback, hModel, X, K, M, S, lamda, idx)
    [feedback, ~, ~] = align_circular_signal(feedback, filter(hModel, 1, x));
    recovered = get_OPLS(K, M, X, x, hModel, feedback, S, lamda);
    [recovered, ~, ~] = align_circular_signal(recovered, target);
    nmseDb = get_NMSE(recovered(idx), target(idx));
end

function writeReport(path, r)
    fid = fopen(path, 'w');
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '# 理想线性区 x,x_f 相位校准\n\n');
    fprintf(fid, '令 x_f=H1(x)，模拟 PA 严格线性且无残余动态相位的校准记录。H1 为 %d 阶椭圆 IIR，85/90 MHz，波纹 0.1 dB；H2 45/48 MHz。使用同一 38400 点实测输入，Fs=1 GHz；只由 x,x_f 选择可靠频点并估相位。该实验并非真实低功率 PA 测量。\n\n', r.h1Order);
    fprintf(fid, 'Welch 4096 点 Hann 窗、2048 点重叠、8192 点 FFT、5 点平滑；门限为输入功率 ≥ 峰值 −30 dB、相干度 ≥ 0.8。实际频率 %.2f 至 %.2f MHz。相位 FIR 4001 抽头、边缘 5 MHz 归零；各组均为纯相位模型与 H2 覆盖。\n\n', ...
        min(r.frequencyMHz), max(r.frequencyMHz));
    fprintf(fid, '| 相位来源 | 对 H1 真值误差 RMS (rad) | 度 |\n|---|---:|---:|\n');
    names = {'x,y_f', 'x,x_f', 'y,y_f 离线算子', 'H1 真值'};
    for k = 1:4
        fprintf(fid, '| %s | %.6f | %.4f |\n', names{k}, ...
            r.phaseErrorRmsRad(k), r.phaseErrorRmsRad(k) * 180 / pi);
    end
    fprintf(fid, '\n将校准相位应用于实测 PA 的 y_f，LD-PLS 参数 K=7、M=5、S=0.3、lamda=1，统一在样本 %d–%d 评分；真实 y 不进入 x,x_f 相位估计。\n\n', r.evaluationIndices);
    fprintf(fid, '| 方法 | 真实 PA 全带恢复 NMSE (dB) |\n|---|---:|\n');
    for k = 1:numel(r.methodLabels)
        fprintf(fid, '| %s | %.2f |\n', strrep(r.methodLabels{k}, '|', '；'), r.nmseDb(k));
    end
    fprintf(fid, '\n线性链路自身恢复 x 的 NMSE：仅 H2 %.2f dB；x,x_f 校准相位 %.2f dB；真实 H1 相位 %.2f dB。\n\n', ...
        r.linearChainNmseDb);
    fprintf(fid, 'x,x_f 校准相对 x,y_f 校准的真实 PA 恢复 NMSE 变化 %.2f dB，与 y,y_f 离线参考相差 %.2f dB。\n\n', ...
        r.nmseDb(3) - r.nmseDb(2), r.nmseDb(3) - r.nmseDb(4));
    fprintf(fid, '该设置完全去除了 PA 本身的相位；真实 PA 即使回退到线性区，也可能保留线性记忆相位，实际低功率校准仍需验证这部分残余。\n');
end
