function results = verify_adaptive_direction_penalty()
% Let the full-band energy penalty change with each extracted PLS direction.
% Select alpha on validation data and compare with the fixed-mu experiment.

fixed = verify_pls_direction_regularization();
projectDir = fileparts(fileparts(mfilename('fullpath')));
dataX = load(fullfile(projectDir, 'data', 'x.mat'), 'x');
dataY = load(fullfile(projectDir, 'data', 'y0.mat'), 'y0');
x = dataX.x(:);
y0 = dataY.y0(:);

K = 7;
M = 5;
lambda = 10;
maxO = round(0.5 * K * (M + 1));
alphaList = [0 logspace(-8, 2, 21)];
ratioPowerList = [0.5 1];
kappaGMax = 10^4.5;
bandwidths = fixed.bandwidths;
nBandwidths = numel(bandwidths);
nAlpha = numel(alphaList);
nPowers = numel(ratioPowerList);
train = 1:floor(0.6 * numel(x));
validation = floor(0.6 * numel(x)) + (1:floor(0.2 * numel(x)));
test = (validation(end) + 1):numel(x);
X = geth_even_1d(x, K, M);

validationNmse = zeros(nBandwidths, maxO, nAlpha, nPowers);
testNmse = zeros(nBandwidths, maxO, nAlpha, nPowers);
deltaByStep = zeros(nBandwidths, maxO, nAlpha, nPowers);
kappaByStep = zeros(nBandwidths, maxO, nAlpha, nPowers);
energyRatioByStep = zeros(nBandwidths, maxO, nAlpha, nPowers);
muByStep = zeros(nBandwidths, maxO, nAlpha, nPowers);
adaptiveO = zeros(nBandwidths, nAlpha, nPowers);
adaptiveValidation = zeros(nBandwidths, nAlpha, nPowers);
adaptiveTest = zeros(nBandwidths, nAlpha, nPowers);

for ib = 1:nBandwidths
    BW = bandwidths(ib);
    h1 = design_filter([BW + 5 BW + 10], 0.5, 1e-4);
    h2 = design_filter([BW BW + 3], 0.1, 1e-6);
    y_f = filter(h1, 1, y0);
    [y_f, ~, ~] = align_circular_signal(y_f, x);
    y_bf = filter(h2, 1, y_f);
    H_x = filter(h2, 1, X);

    xBandLinear = H_x(:, 1);
    wLinear = xBandLinear(train) \ y_bf(train);
    pBand = xBandLinear(train) \ H_x(train, 2:end);
    pFull = X(train, 1) \ X(train, 2:end);
    A = H_x(:, 2:end) - xBandLinear * pBand;
    U = X(:, 2:end) - X(:, 1) * pFull;
    yResidual = y_bf - xBandLinear * wLinear;

    ATrain = A(train, :);
    yTrain = yResidual(train);
    scaleA = std(ATrain);
    scaleY = std(yTrain);
    Z = (ATrain - mean(ATrain)) ./ scaleA;
    zY = (yTrain - mean(yTrain)) ./ scaleY;
    UTrain = (U(train, :) - mean(U(train, :))) ./ scaleA;
    fullGram = UTrain' * UTrain;
    bandGram = Z' * Z;

    for ip = 1:nPowers
        for ia = 1:nAlpha
            [BByO, deltaList, kappaList, ratioList, muList] = ...
                fit_adaptive_penalty_pls(Z, zY, fullGram, bandGram, ...
                scaleA, scaleY, maxO, lambda, alphaList(ia), ratioPowerList(ip));
            yReconstructed = U * BByO + X(:, 1) * wLinear;
            deltaByStep(ib, :, ia, ip) = deltaList;
            kappaByStep(ib, :, ia, ip) = kappaList;
            energyRatioByStep(ib, :, ia, ip) = ratioList;
            muByStep(ib, :, ia, ip) = muList;
            for O = 1:maxO
                validationNmse(ib, O, ia, ip) = get_NMSE( ...
                    yReconstructed(validation, O), y0(validation));
                testNmse(ib, O, ia, ip) = get_NMSE( ...
                    yReconstructed(test, O), y0(test));
            end
            adaptiveO(ib, ia, ip) = select_order(kappaList, maxO, kappaGMax);
            adaptiveValidation(ib, ia, ip) = validationNmse( ...
                ib, adaptiveO(ib, ia, ip), ia, ip);
            adaptiveTest(ib, ia, ip) = testNmse( ...
                ib, adaptiveO(ib, ia, ip), ia, ip);
        end
    end
end

baselineDifference = max(abs( ...
    validationNmse(:, :, 1, :) - fixed.validationNmse(:, :, 1)), [], 'all');
assert(baselineDifference < 1e-8, ...
    'The alpha=0 adaptive implementation does not match original LD-PLS.');

[bestValidation, selectedAlphaIndex] = min(adaptiveValidation, [], 2);
selectedAlpha = zeros(nBandwidths, nPowers);
selectedO = zeros(nBandwidths, nPowers);
selectedTest = zeros(nBandwidths, nPowers);
for ip = 1:nPowers
    for ib = 1:nBandwidths
        ia = selectedAlphaIndex(ib, 1, ip);
        selectedAlpha(ib, ip) = alphaList(ia);
        selectedO(ib, ip) = adaptiveO(ib, ia, ip);
        selectedTest(ib, ip) = adaptiveTest(ib, ia, ip);
    end
end
globalScore = squeeze(mean(adaptiveValidation, 1));
[~, globalIndex] = min(globalScore, [], 'all', 'linear');
[globalAlphaIndex, globalPowerIndex] = ind2sub([nAlpha nPowers], globalIndex);
globalAlpha = alphaList(globalAlphaIndex);
globalPower = ratioPowerList(globalPowerIndex);
globalO = adaptiveO(:, globalAlphaIndex, globalPowerIndex);
globalTest = adaptiveTest(:, globalAlphaIndex, globalPowerIndex);

fprintf('Adaptive mu_k = alpha * (E_full/E_band)^power, stop kappa=%g\n', kappaGMax);
fprintf('BW   fixed mu/O/test      sqrt ratio alpha/O/test   ratio alpha/O/test\n');
for ib = 1:nBandwidths
    fprintf('%3d   %7g / %2d / %8.2f    %8g / %2d / %8.2f    %8g / %2d / %8.2f\n', ...
        bandwidths(ib), fixed.adaptiveSelectedMu(ib, end), ...
        fixed.adaptiveSelectedO(ib, end), fixed.adaptiveSelectedTestNmse(ib, end), ...
        selectedAlpha(ib, 1), selectedO(ib, 1), selectedTest(ib, 1), ...
        selectedAlpha(ib, 2), selectedO(ib, 2), selectedTest(ib, 2));
end
fprintf('Shared alpha=%g, power=%g: O=[%s], test=[%s] dB\n', ...
    globalAlpha, globalPower, num2str(globalO(:).'), num2str(globalTest(:).', ' %.2f'));

results.bandwidths = bandwidths;
results.alphaList = alphaList;
results.ratioPowerList = ratioPowerList;
results.kappaGMax = kappaGMax;
results.validationNmse = validationNmse;
results.testNmse = testNmse;
results.deltaByStep = deltaByStep;
results.kappaByStep = kappaByStep;
results.energyRatioByStep = energyRatioByStep;
results.muByStep = muByStep;
results.adaptiveO = adaptiveO;
results.adaptiveValidation = adaptiveValidation;
results.adaptiveTest = adaptiveTest;
results.selectedAlpha = selectedAlpha;
results.selectedO = selectedO;
results.selectedValidation = bestValidation;
results.selectedTest = selectedTest;
results.globalAlpha = globalAlpha;
results.globalPower = globalPower;
results.globalO = globalO;
results.globalTest = globalTest;
results.fixed = fixed;
end

function [BByO, deltaList, kappaList, ratioList, muList] = ...
    fit_adaptive_penalty_pls(Z, zY, fullGram, bandGram, ...
    scaleA, scaleY, maxO, lambda, alpha, ratioPower)
% Increase the penalty when the unregularized direction is poorly visible.

nFeatures = size(Z, 2);
identity = eye(nFeatures);
W = zeros(nFeatures, maxO);
P = zeros(nFeatures, maxO);
scores = zeros(size(Z, 1), maxO);
scoreGram = zeros(maxO);
deltaList = zeros(1, maxO);
kappaList = zeros(1, maxO);
ratioList = zeros(1, maxO);
muList = zeros(1, maxO);
ZResidual = Z;
yResidual = zY;
previousEnergy = real(yResidual' * yResidual);

for k = 1:maxO
    correlation = ZResidual' * yResidual;
    w0 = correlation / norm(correlation);
    fullEnergy = real(w0' * fullGram * w0);
    bandEnergy = real(w0' * bandGram * w0);
    ratioList(k) = max(fullEnergy, 0) / max(bandEnergy, eps);
    muList(k) = alpha * ratioList(k)^ratioPower;
    if alpha == 0
        w = correlation;
    else
        visibleScale = real(trace(bandGram)) / nFeatures;
        penalty = (fullGram + fullGram') / (2 * visibleScale);
        w = (identity + muList(k) * penalty) \ correlation;
    end
    w = w / norm(w);

    t = ZResidual * w;
    p = ZResidual' * t / (t' * t);
    q = yResidual' * t / (t' * t);
    newYResidual = yResidual - t * q;
    newEnergy = real(newYResidual' * newYResidual);
    deltaList(k) = (previousEnergy - newEnergy) / (previousEnergy + eps);
    scores(:, k) = t;
    scoreGram(1:k-1, k) = scores(:, 1:k-1)' * t;
    scoreGram(k, 1:k-1) = scoreGram(1:k-1, k)';
    scoreGram(k, k) = t' * t;
    kappaList(k) = cond(scoreGram(1:k, 1:k));
    transformation = identity - w * p';
    fullGram = transformation' * fullGram * transformation;
    bandGram = transformation' * bandGram * transformation;
    ZResidual = ZResidual - t * p';
    yResidual = newYResidual;
    previousEnergy = newEnergy;
    W(:, k) = w;
    P(:, k) = p;
end

BByO = zeros(nFeatures, maxO);
for O = 1:maxO
    R = W(:, 1:O) / (P(:, 1:O)' * W(:, 1:O));
    T = Z * R;
    C = (T' * T + lambda * eye(O)) \ (T' * zY);
    BByO(:, O) = (R * C) ./ scaleA' * scaleY;
end
end

function O = select_order(kappaList, maxO, kappaGMax)
firstRejected = find(kappaList(2:maxO) >= kappaGMax, 1);
if isempty(firstRejected)
    O = maxO;
else
    O = firstRejected;
end
end

function h = design_filter(edges, rippleDb, stopRipple)
fs = 1000;
passRipple = (10^(rippleDb / 20) - 1) / (10^(rippleDb / 20) + 1);
[order, frequency, amplitude, weights] = firpmord( ...
    edges / (fs / 2), [1 0], [passRipple stopRipple]);
h = firpm(order, frequency, amplitude, weights, {20});
end
