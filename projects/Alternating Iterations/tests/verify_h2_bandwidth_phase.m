function results = verify_h2_bandwidth_phase(h1Mode)
%VERIFY_H2_BANDWIDTH_PHASE Fixed-H1 bandwidth and phase-control experiment.
% Run from the repository root after startup. H1 truth is used only to
% synthesize observations, diagnose phase, and construct the oracle model.

    if nargin == 0
        h1Mode = 'mixed';
    end
    h1Mode = validatestring(h1Mode, {'mixed', 'phase_only'});

    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    if strcmp(h1Mode, 'phase_only')
        resultFolder = 'h2_phase_only';
    else
        resultFolder = 'h2_bandwidth_phase';
    end
    outDir = fullfile(root, 'projects', 'Alternating Iterations', ...
        'results', resultFolder);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    Fs = 1000;                 % MHz; one sample is one ns
    bandwidths = (45:5:80).'; % positive-frequency passband edges, MHz
    K = 7;
    M = 5;
    S = 0.3;
    lamda = 1;
    operatorR = 250;
    operatorEta = 1e-4;
    Dpass2 = 0.0057501127785;
    Dstop2 = 1e-4;
    guard = 5000;              % common fit-diagnostic and evaluation window
    nImpulse = 3000;

    dataDir = fullfile(root, 'datasets', '测试', '2026.6.10', 'test1');
    d = load(fullfile(dataDir, 'x_v0.6_v1.mat'));
    x = d.x(:);
    d = load(fullfile(dataDir, 'y0_v0.6_v1.mat'));
    y0 = d.y0(:);
    N = numel(x);
    idx = (guard + 1):(N - guard);
    X = geth_even_1d(x, K, M);

    if strcmp(h1Mode, 'phase_only')
        % Two stable real all-pass biquads: exactly unit magnitude, with
        % nonlinear phase strength comparable to the mixed elliptic case.
        r = 0.85;
        theta = 2 * pi * 70 / Fs;
        aStage = [1, -2 * r * cos(theta), r^2];
        bStage = fliplr(aStage);
        a1 = conv(aStage, aStage);
        b1 = conv(bStage, bStage);
        nH1 = numel(a1) - 1;
    else
        [nH1, wnH1] = ellipord(85 / (Fs / 2), 90 / (Fs / 2), 0.5, 80);
        [b1, a1] = ellip(nH1, 0.5, 80, wnH1, 'low');
    end
    assert(max(abs(roots(a1))) < 1, 'H1 is not stable.');
    h1 = impz(b1, a1, nImpulse);
    hTail = impz(b1, a1, 2 * nImpulse);
    impulseTailRatio = norm(hTail(nImpulse + 1:end)) / norm(hTail);
    assert(impulseTailRatio < 1e-3, 'H1 oracle FIR truncation is too large.');

    nFreq = 16384;
    [H1, f] = freqz(b1, a1, nFreq, Fs);
    h1MagnitudeDeviation = max(abs(abs(H1) - 1));
    if strcmp(h1Mode, 'phase_only')
        assert(h1MagnitudeDeviation < 1e-8, 'Phase-only H1 has magnitude error.');
        h1OneDbEdgeMHz = NaN;
    else
        h1OneDbEdgeMHz = f(find(abs(H1) <= 10^(-1 / 20), 1));
    end
    phi = unwrap(angle(H1));
    tauG = grpdelay(b1, a1, nFreq, Fs); % samples = ns at 1 GHz
    omega = 2 * pi * f / Fs;
    h2s = cell(numel(bandwidths), 1);
    phaseRms = zeros(numel(bandwidths), 1);
    groupDelayStdNs = phaseRms;
    h2Order = phaseRms;
    h2GroupDelayNs = phaseRms;
    h2OneDbEdgeMHz = phaseRms;
    h2PassRippleDb = phaseRms;
    h2StopMaxDb = phaseRms;
    singularValues = zeros(numel(bandwidths), size(X, 2));
    effectiveRank = phaseRms;
    conditionNumber = phaseRms;
    cover = initMetrics(numel(bandwidths));
    calibrated = initMetrics(numel(bandwidths));
    oracle = initMetrics(numel(bandwidths));
    operatorFitDb = phaseRms;
    latentCount = round(S * size(X, 2)) * ones(size(bandwidths));

    % One and only one measured PA record and H1 realization for experiment A.
    yH1 = filter(b1, a1, y0);
    [yH1Aligned, xAligned, h1Alignment] = align_circular_signal(yH1, x);
    XAligned = geth_even_1d(xAligned, K, M);
    for j = 1:numel(bandwidths)
        bw = bandwidths(j);
        [ord, fo, ao, weights] = firpmord([bw, bw + 3] / (Fs / 2), ...
            [1 0], [Dpass2 Dstop2]);
        h2 = firpm(ord, fo, ao, weights, {20});
        assert(norm(h2(:) - flipud(h2(:))) < 1e-10 * norm(h2), ...
            'H2 must be symmetric.');
        h2s{j} = h2(:);
        h2Order(j) = ord;
        h2GroupDelayNs(j) = ord / 2;
        H2 = freqz(h2, 1, nFreq, Fs);
        h2OneDbEdgeMHz(j) = f(find(abs(H2) <= 10^(-1 / 20), 1));
        pass = f <= bw;
        stop = f >= bw + 3;
        h2PassRippleDb(j) = 20 * log10(max(abs(H2(pass))) / min(abs(H2(pass))));
        h2StopMaxDb(j) = 20 * log10(max(abs(H2(stop))));

        reliable = pass & abs(H1) >= 10^(-1 / 20);
        fit = [ones(nnz(reliable), 1), omega(reliable)] \ phi(reliable);
        residualPhase = phi(reliable) - [ones(nnz(reliable), 1), omega(reliable)] * fit;
        phaseRms(j) = sqrt(mean(residualPhase.^2));
        groupDelayStdNs(j) = std(tauG(reliable));

        Hx = geth_even_BLMP(xAligned, K, M, h2);
        Hx = Hx(idx, :);
        Hx = Hx ./ vecnorm(Hx);
        sv = svd(Hx, 'econ');
        singularValues(j, :) = sv(:).';
        effectiveRank(j) = nnz(sv >= 1e-3 * sv(1));
        conditionNumber(j) = sv(1) / sv(end);

        yf = filter(h2, 1, yH1Aligned);
        [cover(j), yCover] = fitCover(xAligned, XAligned, y0, yf, h2, ...
            K, M, S, lamda, idx);
        if strcmp(h1Mode, 'mixed')
            [calibrated(j), operatorFitDb(j)] = fitCalibration( ...
                xAligned, XAligned, y0, yf, h2, yCover, ...
                K, M, S, lamda, operatorR, operatorEta, idx);
        else
            operatorFitDb(j) = NaN;
        end

        hOracle = conv(h1(:), h2(:));
        linearReference = filter(hOracle, 1, xAligned);
        [yfOracle, ~, refAlignment] = align_circular_signal(yf, linearReference);
        [yOracle, ~] = get_OPLS(K, M, XAligned, xAligned, ...
            hOracle, yfOracle, S, lamda);
        oracle(j) = scoreMethod(yOracle, y0, yfOracle, hOracle, idx, refAlignment);
        if strcmp(h1Mode, 'phase_only')
            fprintf('H2 %2d MHz: cover %.2f, oracle %.2f dB\n', ...
                bw, cover(j).nmseDb, oracle(j).nmseDb);
        else
            fprintf('H2 %2d MHz: cover %.2f, amplitude %.2f, oracle %.2f dB\n', ...
                bw, cover(j).nmseDb, calibrated(j).nmseDb, oracle(j).nmseDb);
        end
    end

    % Exact offline magnitude/phase separation on a zero-padded DFT grid.
    alphas = [0 0.5 1 1.5];
    phaseControlNmse = zeros(numel(alphas), numel(bandwidths));
    phaseControlRms = phaseControlNmse;
    nDft = 2^nextpow2(3 * N);
    Hwhole = freqz(b1, a1, nDft, 'whole');
    Hpositive = Hwhole(1:nDft / 2 + 1);
    phiPositive = unwrap(angle(Hpositive));
    wPositive = 2 * pi * (0:nDft / 2).' / nDft;
    in80 = wPositive <= 2 * pi * 80 / Fs;
    linearSlope = wPositive(in80) \ phiPositive(in80);
    phiLinear = wPositive * linearSlope;
    yPadded = [zeros(N, 1); y0; zeros(N, 1)];
    for k = 1:numel(alphas)
        alpha = alphas(k);
        phasePositive = phiLinear + alpha * (phiPositive - phiLinear);
        HpositiveControl = abs(Hpositive) .* exp(1j * phasePositive);
        HpositiveControl([1 end]) = real(HpositiveControl([1 end]));
        Hcontrol = [HpositiveControl; conj(flipud(HpositiveControl(2:end - 1)))];
        yControlPadded = ifft(fft(yPadded, nDft) .* Hcontrol);
        yControl = yControlPadded(N + (1:N));
        [yControlAligned, xControlAligned] = align_circular_signal(yControl, x);
        XControl = geth_even_1d(xControlAligned, K, M);
        for j = 1:numel(bandwidths)
            bw = bandwidths(j);
            mask = wPositive <= 2 * pi * bw / Fs;
            v = phasePositive(mask);
            w = wPositive(mask);
            p = [ones(nnz(mask), 1), w] \ v;
            phaseControlRms(k, j) = rms(v - [ones(nnz(mask), 1), w] * p);
            yf = filter(h2s{j}, 1, yControlAligned);
            [m, yCover] = fitCover(xControlAligned, XControl, y0, yf, ...
                h2s{j}, K, M, S, lamda, idx);
            if strcmp(h1Mode, 'mixed')
                [m, ~] = fitCalibration(xControlAligned, XControl, y0, yf, ...
                    h2s{j}, yCover, K, M, S, lamda, ...
                    operatorR, operatorEta, idx);
            end
            phaseControlNmse(k, j) = m.nmseDb;
        end
    end

    % Stable second-order all-pass cascade: its magnitude is one by design.
    r = 0.85;
    theta = 2 * pi * 70 / Fs;
    aAp = [1, -2 * r * cos(theta), r^2];
    bAp = fliplr(aAp);
    assert(max(abs(roots(aAp))) < 1);
    Hap = freqz(bAp, aAp, nFreq, Fs);
    assert(max(abs(abs(Hap) - 1)) < 1e-8);
    allpassBands = [50 65 80];
    if strcmp(h1Mode, 'phase_only')
        stageCounts = 0:3;
    else
        stageCounts = 0:2;
    end
    allpassNmse = zeros(numel(stageCounts), numel(allpassBands));
    Hmagnitude = [abs(Hpositive); flipud(abs(Hpositive(2:end - 1)))];
    yMagnitudePadded = ifft(fft(yPadded, nDft) .* Hmagnitude);
    yAllpass = yMagnitudePadded(N + (1:N));
    for stage = stageCounts
        if stage > 0
            yAllpass = filter(bAp, aAp, yAllpass);
        end
        [yApAligned, xApAligned] = align_circular_signal(yAllpass, x);
        XAp = geth_even_1d(xApAligned, K, M);
        for j = 1:numel(allpassBands)
            h2 = h2s{bandwidths == allpassBands(j)};
            yf = filter(h2, 1, yApAligned);
            [m, yCover] = fitCover(xApAligned, XAp, y0, yf, ...
                h2, K, M, S, lamda, idx);
            if strcmp(h1Mode, 'mixed')
                [m, ~] = fitCalibration(xApAligned, XAp, y0, yf, ...
                    h2, yCover, K, M, S, lamda, operatorR, operatorEta, idx);
            end
            allpassNmse(stage + 1, j) = m.nmseDb;
        end
    end

    if strcmp(h1Mode, 'phase_only')
        summary = table(bandwidths, h2Order, h2GroupDelayNs, h2OneDbEdgeMHz, ...
            h2PassRippleDb, h2StopMaxDb, phaseRms, groupDelayStdNs, ...
            [cover.nmseDb].', [oracle.nmseDb].', ...
            [cover.fitDb].', [oracle.fitDb].', ...
            effectiveRank, conditionNumber, latentCount, ...
            'VariableNames', {'BandwidthMHz', 'H2Order', 'H2DelayNs', ...
            'H2OneDbEdgeMHz', 'H2PassRippleDb', 'H2StopMaxDb', ...
            'PhaseRmsRad', 'GroupDelayStdNs', 'CoverNmseDb', ...
            'OracleNmseDb', 'CoverFitDb', 'OracleFitDb', ...
            'EffectiveRank', 'ConditionNumber', 'LatentCount'});
        detailed = table(bandwidths, ...
            [cover.nmseDb].', [cover.rawNmseDb].', ...
            [cover.optimalGainNmseDb].', [oracle.nmseDb].', ...
            [oracle.rawNmseDb].', [oracle.optimalGainNmseDb].', ...
            'VariableNames', {'BandwidthMHz', 'CoverNormalizedDb', ...
            'CoverRawDb', 'CoverOptimalGainDb', 'OracleNormalizedDb', ...
            'OracleRawDb', 'OracleOptimalGainDb'});
    else
        summary = table(bandwidths, h2Order, h2GroupDelayNs, h2OneDbEdgeMHz, ...
            h2PassRippleDb, h2StopMaxDb, phaseRms, groupDelayStdNs, ...
            [cover.nmseDb].', [calibrated.nmseDb].', [oracle.nmseDb].', ...
            [cover.fitDb].', [calibrated.fitDb].', [oracle.fitDb].', ...
            operatorFitDb, effectiveRank, conditionNumber, latentCount, ...
            'VariableNames', {'BandwidthMHz', 'H2Order', 'H2DelayNs', ...
            'H2OneDbEdgeMHz', 'H2PassRippleDb', 'H2StopMaxDb', ...
            'PhaseRmsRad', 'GroupDelayStdNs', 'CoverNmseDb', ...
            'AmplitudeNmseDb', 'OracleNmseDb', 'CoverFitDb', ...
            'AmplitudeFitDb', 'OracleFitDb', 'OperatorFitDb', ...
            'EffectiveRank', 'ConditionNumber', 'LatentCount'});
        detailed = table(bandwidths, ...
            [cover.nmseDb].', [cover.rawNmseDb].', [cover.optimalGainNmseDb].', ...
            [calibrated.nmseDb].', [calibrated.rawNmseDb].', ...
            [calibrated.optimalGainNmseDb].', [oracle.nmseDb].', ...
            [oracle.rawNmseDb].', [oracle.optimalGainNmseDb].', ...
            'VariableNames', {'BandwidthMHz', 'CoverNormalizedDb', ...
            'CoverRawDb', 'CoverOptimalGainDb', 'AmplitudeNormalizedDb', ...
            'AmplitudeRawDb', 'AmplitudeOptimalGainDb', ...
            'OracleNormalizedDb', 'OracleRawDb', 'OracleOptimalGainDb'});
    end
    writetable(summary, fullfile(outDir, 'summary.csv'));
    writetable(detailed, fullfile(outDir, 'nmse_details.csv'));
    controlTable = array2table(phaseControlNmse.', ...
        'VariableNames', {'Phase0', 'Phase05', 'Phase1', 'Phase15'});
    controlTable = addvars(controlTable, bandwidths, 'Before', 1, ...
        'NewVariableNames', 'BandwidthMHz');
    writetable(controlTable, fullfile(outDir, 'phase_controls.csv'));
    phaseTable = array2table(phaseControlRms.', ...
        'VariableNames', {'Phase0', 'Phase05', 'Phase1', 'Phase15'});
    phaseTable = addvars(phaseTable, bandwidths, 'Before', 1, ...
        'NewVariableNames', 'BandwidthMHz');
    writetable(phaseTable, fullfile(outDir, 'phase_control_rms.csv'));

    results = struct('summary', summary, 'cover', cover, ...
        'calibrated', calibrated, 'oracle', oracle, ...
        'singularValues', singularValues, 'h1Order', nH1, ...
        'h1Numerator', b1, 'h1Denominator', a1, ...
        'h1FrequencyMHz', f, 'h1Response', H1, ...
        'h1UnwrappedPhase', phi, 'h1GroupDelayNs', tauG, ...
        'h1Mode', h1Mode, 'h1MagnitudeDeviation', h1MagnitudeDeviation, ...
        'h1OneDbEdgeMHz', h1OneDbEdgeMHz, ...
        'h2Coefficients', {h2s}, ...
        'h1ImpulseTailRatio', impulseTailRatio, ...
        'h1Alignment', h1Alignment, 'phaseStrength', alphas, ...
        'phaseControlNmse', phaseControlNmse, ...
        'phaseControlRms', phaseControlRms, ...
        'allpassBands', allpassBands, 'allpassNmse', allpassNmse, ...
        'allpassStageCounts', stageCounts, ...
        'evaluationIndices', [idx(1) idx(end)], ...
        'dataFile', dataDir);
    if strcmp(h1Mode, 'phase_only')
        results = rmfield(results, 'calibrated');
    end
    save(fullfile(outDir, 'results.mat'), 'results');
    makeFigures(outDir, f, H1, phi, tauG, bandwidths, ...
        phaseRms, groupDelayStdNs, cover, calibrated, oracle, ...
        singularValues, effectiveRank, conditionNumber, ...
        alphas, phaseControlNmse, allpassBands, allpassNmse, ...
        stageCounts, h1Mode);
    disp(summary);
end

function m = initMetrics(n)
    m = repmat(struct('nmseDb', NaN, 'rawNmseDb', NaN, ...
        'optimalGainNmseDb', NaN, 'fitDb', NaN, 'alignment', NaN), n, 1);
end

function [m, yCover] = fitCover(x, X, y0, yf, h2, K, M, S, lamda, idx)
    [yCover, ~] = get_OPLS(K, M, X, x, h2, yf, S, lamda);
    m = scoreMethod(yCover, y0, yf, h2, idx, NaN);
end

function [m, operatorFitDb] = fitCalibration( ...
        x, X, y0, yf, h2, yCover, K, M, S, lamda, R, eta, idx)
    [hModel, ~, ~, operatorFitDb] = estimate_amplitude_operator( ...
        yCover, yf, h2, R, eta);
    linearReference = filter(hModel, 1, x);
    [yfAligned, ~, feedbackAlignment] = align_circular_signal(yf, linearReference);
    [yCal, ~] = get_OPLS(K, M, X, x, hModel, yfAligned, S, lamda);
    m = scoreMethod(yCal, y0, yfAligned, hModel, idx, feedbackAlignment);
end

function m = scoreMethod(yHat, y0, yf, hModel, idx, feedbackAlignment)
    [yHatAligned, ~, evaluationAlignment] = align_circular_signal(yHat, y0);
    reference = y0(idx);
    estimate = yHatAligned(idx);
    m.nmseDb = get_NMSE(estimate, reference);
    m.rawNmseDb = 10 * log10(norm(estimate - reference)^2 / norm(reference)^2);
    c = (estimate' * reference) / (estimate' * estimate);
    m.optimalGainNmseDb = 10 * log10(norm(c * estimate - reference)^2 / norm(reference)^2);
    modelFeedback = filter(hModel, 1, yHat);
    m.fitDb = 10 * log10(norm(modelFeedback(idx) - yf(idx))^2 / norm(yf(idx))^2);
    m.alignment = [feedbackAlignment, evaluationAlignment];
end

function makeFigures(outDir, f, H1, phi, tauG, bw, ePhi, tauStd, ...
        cover, cal, oracle, sv, rankEff, condNum, alphas, nmseControl, ...
        allpassBands, allpassNmse, stageCounts, h1Mode)
    fig = figure('Visible', 'off');
    tiledlayout(3, 1);
    nexttile; plot(f, 20 * log10(abs(H1))); ylabel('|H1| (dB)'); ylim([-100 5]);
    addEdges(bw); title('H1 and H2 passband edges'); grid on;
    nexttile; plot(f, phi); ylabel('Unwrapped phase (rad)');
    addEdges(bw); grid on;
    nexttile; plot(f, tauG); ylabel('Group delay (ns)'); xlabel('Frequency (MHz)');
    addEdges(bw); xlim([0 100]); grid on;
    exportgraphics(fig, fullfile(outDir, 'figure1_h1_response.png')); close(fig);

    fig = figure('Visible', 'off');
    yyaxis left; plot(bw, ePhi, '-o'); ylabel('Residual phase RMS (rad)');
    yyaxis right; plot(bw, tauStd, '-s'); ylabel('Group delay std (ns)');
    xlabel('H2 passband edge (MHz)'); grid on;
    exportgraphics(fig, fullfile(outDir, 'figure2_phase_vs_bandwidth.png')); close(fig);

    fig = figure('Visible', 'off');
    if strcmp(h1Mode, 'phase_only')
        plot(bw, [cover.nmseDb], '-o', bw, [oracle.nmseDb], '-^');
        legend('Cover', 'Known complex H1H2', 'Location', 'best');
    else
        plot(bw, [cover.nmseDb], '-o', bw, [cal.nmseDb], '-s', ...
            bw, [oracle.nmseDb], '-^');
        legend('Cover', 'One amplitude calibration', 'Known complex H1H2', ...
            'Location', 'best');
    end
    xlabel('H2 passband edge (MHz)'); ylabel('Full-band NMSE (dB)');
    grid on;
    exportgraphics(fig, fullfile(outDir, 'figure3_fullband_nmse.png')); close(fig);

    fig = figure('Visible', 'off');
    tiledlayout(2, 1);
    nexttile; semilogy(bw, condNum, '-o'); ylabel('Column-normalized condition number'); grid on;
    nexttile; plot(bw, rankEff, '-s'); ylabel('Effective rank');
    xlabel('H2 passband edge (MHz)'); grid on;
    exportgraphics(fig, fullfile(outDir, 'figure4_matrix_condition.png')); close(fig);
    fig = figure('Visible', 'off');
    semilogy(1:size(sv, 2), sv ./ sv(:, 1), 'LineWidth', 1);
    xlabel('Singular value index'); ylabel('Relative singular value'); grid on;
    exportgraphics(fig, fullfile(outDir, 'singular_value_spectra.png')); close(fig);

    fig = figure('Visible', 'off');
    plot(bw, nmseControl.', '-o');
    xlabel('H2 passband edge (MHz)'); ylabel('Full-band NMSE (dB)');
    legend(compose('Phase strength %.1f', alphas), 'Location', 'best'); grid on;
    exportgraphics(fig, fullfile(outDir, 'figure5_phase_strength.png')); close(fig);
    fig = figure('Visible', 'off');
    plot(allpassBands, allpassNmse.', '-o');
    xlabel('H2 passband edge (MHz)'); ylabel('Full-band NMSE (dB)');
    legend(compose('%d all-pass stages', stageCounts), ...
        'Location', 'best'); grid on;
    exportgraphics(fig, fullfile(outDir, 'allpass_time_domain_check.png')); close(fig);
end

function addEdges(bw)
    hold on;
    for j = 1:numel(bw)
        xline(bw(j), ':', 'Color', [0.7 0.7 0.7]);
    end
    xlim([0 100]);
end
