function r = diag_high_ripple_adaptive()
%DIAG_HIGH_RIPPLE_ADAPTIVE 分解 1.9 dB 纹波点的幅度估计与重辨识误差。

    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    dataDir = fullfile(root, 'datasets', '测试', '2026.6.10', 'test1');
    d = load(fullfile(dataDir, 'x_v0.6_v1.mat'));
    x = d.x(:);
    d = load(fullfile(dataDir, 'y0_v0.6_v1.mat'));
    y0 = d.y0(:);
    Fs = 1000;
    BW = 80;
    K = 7;
    M = 5;
    lambda = 10;
    T_delta = 0;
    kappaGMax = 10^4.5;

    flatRipple = (10^(0.01 / 20) - 1) / (10^(0.01 / 20) + 1);
    [n1, fo, ao, w] = firpmord([BW + 5, BW + 10] / (Fs / 2), ...
        [1 0], [flatRipple, 1e-4]);
    hFlat = firpm(n1, fo, ao, w, {20});
    tau = round(5 * Fs / (BW + 5));
    epsilon = (10^(1.9 / 20) - 1) / (10^(1.9 / 20) + 1);
    hRipple = [epsilon / 2, zeros(1, tau - 1), 1, ...
        zeros(1, tau - 1), epsilon / 2];
    h1 = conv(hFlat, hRipple);
    [n2, fo, ao, w] = firpmord([BW, BW + 3] / (Fs / 2), ...
        [1 0], [0.0057501127785, 1e-4]);
    h2 = firpm(n2, fo, ao, w, {20});

    yf = filter(h1, 1, y0);
    [yf, xa, ~] = align_circular_signal(yf, x);
    ybf = filter(h2, 1, yf);
    X = geth_even_1d(xa, K, M);
    [yInitial, ~, oInitial] = get_adaptive_LD_PLS( ...
        K, M, X, xa, h2, ybf, lambda, T_delta, kappaGMax);
    hRef = conv(h1, h2);
    ybfRef = align_circular_signal(ybf, filter(hRef, 1, xa));
    [yRef, ~, oRef] = get_adaptive_LD_PLS( ...
        K, M, X, xa, hRef, ybfRef, lambda, T_delta, kappaGMax);
    yRef = align_circular_signal(yRef, y0);

    yCurrent = yInitial;
    iterationNmse = zeros(4, 1);
    iterationAlignedNmse = zeros(4, 1);
    iterationNmse(1) = get_NMSE(yCurrent, y0);
    iterationAlignedNmse(1) = get_NMSE(align_circular_signal(yCurrent, y0), y0);
    operatorFitDb = zeros(3, 1);
    operatorErrorDb = zeros(3, 1);
    orders = zeros(3, 1);
    for i = 1:3
        [hModel, ~, ~, operatorFitDb(i)] = ...
            estimate_amplitude_operator(yCurrent, ybf, h2, 250, 1e-4);
        operatorErrorDb(i) = passband_amplitude_error(hModel, h1, h2, Fs, BW);
        ybfAligned = align_circular_signal(ybf, filter(hModel, 1, xa));
        [yCurrent, ~, orders(i)] = get_adaptive_LD_PLS( ...
            K, M, X, xa, hModel, ybfAligned, lambda, T_delta, kappaGMax);
        iterationNmse(i + 1) = get_NMSE(yCurrent, y0);
        iterationAlignedNmse(i + 1) = get_NMSE( ...
            align_circular_signal(yCurrent, y0), y0);
    end

    [hModelOracle, ~, ~, oracleFitDb] = ...
        estimate_amplitude_operator(y0, ybf, h2, 250, 1e-4);
    oracleOperatorErrorDb = passband_amplitude_error( ...
        hModelOracle, h1, h2, Fs, BW);
    ybfOracle = align_circular_signal(ybf, filter(hModelOracle, 1, xa));
    [yOracle, ~, oOracle] = get_adaptive_LD_PLS( ...
        K, M, X, xa, hModelOracle, ybfOracle, lambda, T_delta, kappaGMax);

    operatorOrders = [60, 100, 250];
    regularizations = [1e-6, 1e-4, 1e-2];
    sweep = zeros(numel(operatorOrders) * numel(regularizations), 5);
    sweepRow = 0;
    for i = 1:numel(operatorOrders)
        for j = 1:numel(regularizations)
            R = operatorOrders(i);
            eta = regularizations(j);
            [hCandidate, ~, ~, fitDb] = ...
                estimate_amplitude_operator(yInitial, ybf, h2, R, eta);
            ybfCandidate = align_circular_signal( ...
                ybf, filter(hCandidate, 1, xa));
            [yCandidate, ~, selectedOrder] = get_adaptive_LD_PLS( ...
                K, M, X, xa, hCandidate, ybfCandidate, ...
                lambda, T_delta, kappaGMax);
            sweepRow = sweepRow + 1;
            sweep(sweepRow, :) = [R, eta, fitDb, ...
                passband_amplitude_error(hCandidate, h1, h2, Fs, BW), ...
                get_NMSE(yCandidate, y0)];
            fprintf('R=%3d eta=%g O=%d fit=%7.2f dB H1err=%5.3f dB NMSE=%7.2f dB\n', ...
                R, eta, selectedOrder, sweep(sweepRow, 3), ...
                sweep(sweepRow, 4), sweep(sweepRow, 5));
        end
    end

    lambdaValues = [1, 3, 10, 30, 100];
    lambdaSweep = zeros(numel(lambdaValues), 6);
    for i = 1:numel(lambdaValues)
        lambdaCandidate = lambdaValues(i);
        [yStart, ~, oStart] = get_adaptive_LD_PLS( ...
            K, M, X, xa, h2, ybf, lambdaCandidate, T_delta, kappaGMax);
        hCandidate = estimate_amplitude_operator(yStart, ybf, h2, 250, 1e-4);
        ybfCandidate = align_circular_signal( ...
            ybf, filter(hCandidate, 1, xa));
        [yCalibrated, ~, oCalibrated] = get_adaptive_LD_PLS( ...
            K, M, X, xa, hCandidate, ybfCandidate, ...
            lambdaCandidate, T_delta, kappaGMax);
        yOracleCandidate = get_adaptive_LD_PLS( ...
            K, M, X, xa, hRef, ybfRef, ...
            lambdaCandidate, T_delta, kappaGMax);
        yOracleCandidate = align_circular_signal(yOracleCandidate, y0);
        lambdaSweep(i, :) = [lambdaCandidate, oStart, ...
            get_NMSE(yStart, y0), oCalibrated, ...
            get_NMSE(yCalibrated, y0), get_NMSE(yOracleCandidate, y0)];
        fprintf('lambda=%3g: Initial %7.2f(O=%d), Iteration %7.2f(O=%d), Ref %7.2f dB\n', ...
            lambdaCandidate, lambdaSweep(i, 3), oStart, ...
            lambdaSweep(i, 5), oCalibrated, lambdaSweep(i, 6));
    end

    r = struct('initialOrder', oInitial, 'referenceOrder', oRef, ...
        'referenceNmse', get_NMSE(yRef, y0), ...
        'iterationNmse', iterationNmse, 'iterationOrders', orders, ...
        'iterationAlignedNmse', iterationAlignedNmse, ...
        'operatorFitDb', operatorFitDb, ...
        'operatorErrorDb', operatorErrorDb, ...
        'oracleOperatorFitDb', oracleFitDb, ...
        'oracleOperatorErrorDb', oracleOperatorErrorDb, ...
        'oracleOperatorOrder', oOracle, ...
        'oracleOperatorNmse', get_NMSE(yOracle, y0), ...
        'parameterSweep', sweep, 'lambdaSweep', lambdaSweep);
    disp(r);
end

function errorDb = passband_amplitude_error(hModel, h1, h2, Fs, BW)
    [hm, f] = freqz(hModel, 1, 8192, Fs);
    h1Response = freqz(h1, 1, 8192, Fs);
    h2Response = freqz(h2, 1, 8192, Fs);
    idx = f <= BW;
    estimated = abs(hm(idx) ./ h2Response(idx));
    actual = abs(h1Response(idx));
    estimated = estimated / estimated(1);
    actual = actual / actual(1);
    errorDb = sqrt(mean((20 * log10(estimated) - 20 * log10(actual)).^2));
end
