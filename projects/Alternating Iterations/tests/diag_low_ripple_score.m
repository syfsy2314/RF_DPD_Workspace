function r = diag_low_ripple_score()
%DIAG_LOW_RIPPLE_SCORE 核对旧版固定潜变量实验的 0.1 dB 纹波点。

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
    S = 0.3;
    lambda = 1;

    flatRipple = (10^(0.01 / 20) - 1) / (10^(0.01 / 20) + 1);
    [n1, fo, ao, w] = firpmord([BW + 5, BW + 10] / (Fs / 2), ...
        [1 0], [flatRipple, 1e-4]);
    hFlat = firpm(n1, fo, ao, w, {20});
    tau = round(5 * Fs / (BW + 5));
    epsilon = (10^(0.1 / 20) - 1) / (10^(0.1 / 20) + 1);
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
    yInitial = get_LD_PLS(K, M, X, xa, h2, ybf, S, lambda);
    [yInitialAligned, ~, ~] = align_circular_signal(yInitial, y0);
    hModel = estimate_amplitude_operator(yInitial, ybf, h2, 250, 1e-4);
    ybfIteration = align_circular_signal(ybf, filter(hModel, 1, xa));
    yIteration = get_LD_PLS(K, M, X, xa, hModel, ybfIteration, S, lambda);
    [yIterationAligned, ~, ~] = align_circular_signal(yIteration, y0);

    hRef = conv(h1, h2);
    [ybfRef, ~, ~] = align_circular_signal(ybf, filter(hRef, 1, xa));
    yRef = get_LD_PLS(K, M, X, xa, hRef, ybfRef, S, lambda);
    [yRef, ~, ~] = align_circular_signal(yRef, y0);
    [yTarget, ~, ~] = align_circular_signal(y0, xa);
    yDirect = X * (X \ yTarget);
    [yDirect, ~, ~] = align_circular_signal(yDirect, y0);

    latentCounts = [13, 20, 30];
    regularizationStrengths = [0, 1];
    nmseGrid = zeros(2, numel(regularizationStrengths), numel(latentCounts));
    modelFilters = {h2, hRef};
    modelFeedback = {ybf, ybfRef};
    for modelIdx = 1:2
        Hx = geth_even_BLMP(xa, K, M, modelFilters{modelIdx});
        for muIdx = 1:numel(regularizationStrengths)
            for latentIdx = 1:numel(latentCounts)
                coefficients = get_energy_regularized_LD_PLS( ...
                    X, Hx, modelFeedback{modelIdx}, lambda, ...
                    regularizationStrengths(muIdx), 0, inf, ...
                    latentCounts(latentIdx));
                yCandidate = X * coefficients;
                [yCandidate, ~, ~] = align_circular_signal(yCandidate, y0);
                nmseGrid(modelIdx, muIdx, latentIdx) = get_NMSE(yCandidate, y0);
            end
        end
    end

    r = struct('initialRaw', get_NMSE(yInitial, y0), ...
        'initialAligned', get_NMSE(yInitialAligned, y0), ...
        'iterationRaw', get_NMSE(yIteration, y0), ...
        'iterationAligned', get_NMSE(yIterationAligned, y0), ...
        'reference', get_NMSE(yRef, y0), ...
        'directFullbandLs', get_NMSE(yDirect, y0), ...
        'latentCounts', latentCounts, ...
        'regularizationStrengths', regularizationStrengths, ...
        'nmseGrid', nmseGrid);
    disp(r);
end
