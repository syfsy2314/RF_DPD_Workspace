function results = verify_pls_direction_regularization()
% Compare original LD-PLS and full-band energy direction regularization.
% Select O and mu on the validation segment; use the test segment only for evaluation.

projectDir = fileparts(fileparts(mfilename('fullpath')));
dataX = load(fullfile(projectDir, 'data', 'x.mat'), 'x');
dataY = load(fullfile(projectDir, 'data', 'y0.mat'), 'y0');
x = dataX.x(:);
y0 = dataY.y0(:);

K = 7;
M = 5;
OList = 1:round(0.5 * K * (M + 1));
referenceO = round(0.3 * K * (M + 1));
checkOList = unique([1 referenceO OList(end)]);
lambda = 10;
fs = 1000;
bandwidths = [80 50 30];
muList = [0 1e-4 1e-3 1e-2 0.1 1 10 30 100 300];
T_delta = 0;
% Threshold selected for energy-regularized directions; see the adaptive sweep
% for the separately calibrated original LD-PLS threshold (1e4).
kappaGMax = 10^4.5;
adaptiveCaps = [referenceO OList(end)];

nTrain = floor(0.6 * numel(x));
nValidation = floor(0.2 * numel(x));
train = 1:nTrain;
validation = nTrain + (1:nValidation);
test = (nTrain + nValidation + 1):numel(x);

X = geth_even_1d(x, K, M);
nBandwidths = numel(bandwidths);
nO = numel(OList);
nMu = numel(muList);
validationNmse = zeros(nBandwidths, nO, nMu);
testNmse = zeros(nBandwidths, nO, nMu);
deltaByStep = zeros(nBandwidths, nO, nMu);
kappaByStep = zeros(nBandwidths, nO, nMu);
baselineDifference = zeros(nBandwidths, 1);
originalO = zeros(nBandwidths, 1);
originalTestNmse = zeros(nBandwidths, 1);
fixedOriginalTestNmse = zeros(nBandwidths, 1);
fixedEnergyMu = zeros(nBandwidths, 1);
fixedEnergyTestNmse = zeros(nBandwidths, 1);
selectedO = zeros(nBandwidths, 1);
selectedMu = zeros(nBandwidths, 1);
selectedTestNmse = zeros(nBandwidths, 1);
adaptiveO = zeros(nBandwidths, numel(adaptiveCaps), nMu);
adaptiveOriginalTestNmse = zeros(nBandwidths, numel(adaptiveCaps));
adaptiveSelectedO = zeros(nBandwidths, numel(adaptiveCaps));
adaptiveSelectedMu = zeros(nBandwidths, numel(adaptiveCaps));
adaptiveSelectedTestNmse = zeros(nBandwidths, numel(adaptiveCaps));
dompBicTestNmse = zeros(nBandwidths, 1);
dompSelectedTestNmse = zeros(nBandwidths, 1);
dompBicOrder = zeros(nBandwidths, 1);
dompSelectedOrder = zeros(nBandwidths, 1);

for ib = 1:nBandwidths
    BW = bandwidths(ib);
    h1 = design_filter([BW + 5 BW + 10], 0.5, 1e-4, fs);
    h2 = design_filter([BW BW + 3], 0.1, 1e-6, fs);

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

    for im = 1:nMu
        [BByO, deltaList, kappaList] = fit_direction_regularized_pls( ...
            Z, zY, fullGram, bandGram, scaleA, scaleY, OList(end), lambda, muList(im));
        yReconstructed = U * BByO + X(:, 1) * wLinear;
        deltaByStep(ib, :, im) = deltaList;
        kappaByStep(ib, :, im) = kappaList;

        for ic = 1:numel(adaptiveCaps)
            adaptiveO(ib, ic, im) = select_adaptive_order( ...
                deltaList, kappaList, adaptiveCaps(ic), T_delta, kappaGMax);
        end

        for io = 1:nO
            validationNmse(ib, io, im) = get_NMSE( ...
                yReconstructed(validation, io), y0(validation));
            testNmse(ib, io, im) = get_NMSE( ...
                yReconstructed(test, io), y0(test));
        end

        if im == 1
            for checkO = checkOList
                reference = get_LD_PLS(K, M, X(train, :), x(train), ...
                    h2, y_bf(train), checkO / (K * (M + 1)), lambda);
                difference = norm( ...
                    yReconstructed(train, checkO) - reference) / norm(reference);
                baselineDifference(ib) = max(baselineDifference(ib), difference);
            end
            assert(baselineDifference(ib) < 1e-6, ...
                'The mu=0 result does not reproduce the original LD-PLS fit.');
        end
    end

    [~, originalIndex] = min(validationNmse(ib, :, 1));
    originalO(ib) = OList(originalIndex);
    originalTestNmse(ib) = testNmse(ib, originalIndex, 1);
    fixedOriginalTestNmse(ib) = testNmse(ib, referenceO, 1);
    [~, fixedMuIndex] = min(reshape(validationNmse(ib, referenceO, :), 1, nMu));
    fixedEnergyMu(ib) = muList(fixedMuIndex);
    fixedEnergyTestNmse(ib) = testNmse(ib, referenceO, fixedMuIndex);
    [~, selectedIndex] = min(reshape(validationNmse(ib, :, :), nO, nMu), [], 'all', 'linear');
    [io, im] = ind2sub([nO, nMu], selectedIndex);
    selectedO(ib) = OList(io);
    selectedMu(ib) = muList(im);
    selectedTestNmse(ib) = testNmse(ib, io, im);

    for ic = 1:numel(adaptiveCaps)
        adaptiveOriginalTestNmse(ib, ic) = testNmse(ib, adaptiveO(ib, ic, 1), 1);
        adaptiveValidationNmse = zeros(1, nMu);
        for im = 1:nMu
            adaptiveValidationNmse(im) = validationNmse(ib, adaptiveO(ib, ic, im), im);
        end
        [~, adaptiveMuIndex] = min(adaptiveValidationNmse);
        adaptiveSelectedO(ib, ic) = adaptiveO(ib, ic, adaptiveMuIndex);
        adaptiveSelectedMu(ib, ic) = muList(adaptiveMuIndex);
        adaptiveSelectedTestNmse(ib, ic) = testNmse( ...
            ib, adaptiveSelectedO(ib, ic), adaptiveMuIndex);
    end

    labels = cell(1, size(H_x, 2));
    [hBic, ~, dompBicOrder(ib), hHistory] = omp_domp( ...
        H_x(train, :), y_bf(train), labels, 'DOMP');
    dompBicTestNmse(ib) = get_NMSE(X(test, :) * hBic, y0(test));
    dompValidationNmse = zeros(1, size(hHistory, 2));
    for n = 1:size(hHistory, 2)
        dompValidationNmse(n) = get_NMSE( ...
            X(validation, :) * hHistory(:, n), y0(validation));
    end
    [~, dompSelectedOrder(ib)] = min(dompValidationNmse);
    dompSelectedTestNmse(ib) = get_NMSE( ...
        X(test, :) * hHistory(:, dompSelectedOrder(ib)), y0(test));
end

fprintf('First rejected energy candidate: BW   O   Delta   kappa\n');
for ib = 1:nBandwidths
    im = find(muList == adaptiveSelectedMu(ib, end), 1);
    rejectedO = adaptiveSelectedO(ib, end) + 1;
    if rejectedO <= OList(end)
        fprintf('%3d   %2d   %8.3g   %8.3g\n', bandwidths(ib), rejectedO, ...
            deltaByStep(ib, rejectedO, im), kappaByStep(ib, rejectedO, im));
    end
end

fprintf('Adaptive stop: T_delta=%g, kappa_G_max=%g\n', T_delta, kappaGMax);
fprintf('BW   cap   original O/test   energy O/mu/test\n');
for ib = 1:nBandwidths
    for ic = 1:numel(adaptiveCaps)
        fprintf('%3d   %2d      %2d / %8.2f       %2d / %5g / %8.2f\n', ...
            bandwidths(ib), adaptiveCaps(ic), adaptiveO(ib, ic, 1), ...
            adaptiveOriginalTestNmse(ib, ic), adaptiveSelectedO(ib, ic), ...
            adaptiveSelectedMu(ib, ic), adaptiveSelectedTestNmse(ib, ic));
    end
end

fprintf('BW   original scan O/test   fixed O=13 original/energy(mu)   energy scan O/mu/test   DOMP order/test\n');
for ib = 1:nBandwidths
    fprintf('%3d   %2d / %8.2f         %8.2f / %8.2f (%5g)        %2d / %5g / %8.2f    %2d / %8.2f\n', ...
        bandwidths(ib), originalO(ib), originalTestNmse(ib), ...
        fixedOriginalTestNmse(ib), fixedEnergyTestNmse(ib), fixedEnergyMu(ib), ...
        selectedO(ib), selectedMu(ib), selectedTestNmse(ib), ...
        dompSelectedOrder(ib), dompSelectedTestNmse(ib));
end

results.bandwidths = bandwidths;
results.OList = OList;
results.muList = muList;
results.validationNmse = validationNmse;
results.testNmse = testNmse;
results.deltaByStep = deltaByStep;
results.kappaByStep = kappaByStep;
results.originalO = originalO;
results.originalTestNmse = originalTestNmse;
results.fixedO = referenceO;
results.fixedOriginalTestNmse = fixedOriginalTestNmse;
results.fixedEnergyMu = fixedEnergyMu;
results.fixedEnergyTestNmse = fixedEnergyTestNmse;
results.selectedO = selectedO;
results.selectedMu = selectedMu;
results.selectedTestNmse = selectedTestNmse;
results.T_delta = T_delta;
results.kappaGMax = kappaGMax;
results.adaptiveCaps = adaptiveCaps;
results.adaptiveO = adaptiveO;
results.adaptiveOriginalTestNmse = adaptiveOriginalTestNmse;
results.adaptiveSelectedO = adaptiveSelectedO;
results.adaptiveSelectedMu = adaptiveSelectedMu;
results.adaptiveSelectedTestNmse = adaptiveSelectedTestNmse;
results.baselineDifference = baselineDifference;
results.dompBicTestNmse = dompBicTestNmse;
results.dompSelectedTestNmse = dompSelectedTestNmse;
results.dompBicOrder = dompBicOrder;
results.dompSelectedOrder = dompSelectedOrder;
end

function h = design_filter(edges, rippleDb, stopRipple, fs)
passRipple = (10^(rippleDb / 20) - 1) / (10^(rippleDb / 20) + 1);
[order, frequency, amplitude, weights] = firpmord( ...
    edges / (fs / 2), [1 0], [passRipple stopRipple]);
h = firpm(order, frequency, amplitude, weights, {20});
end

function [BByO, deltaList, kappaList] = fit_direction_regularized_pls( ...
    Z, zY, fullGram, bandGram, scaleA, scaleY, OMax, lambda, mu)
% Penalize full-band energy in each PLS extraction direction; mu=0 follows pls_B.

nFeatures = size(Z, 2);
W = zeros(nFeatures, OMax);
P = zeros(nFeatures, OMax);
scores = zeros(size(Z, 1), OMax);
scoreGram = zeros(OMax);
deltaList = zeros(1, OMax);
kappaList = zeros(1, OMax);
ZResidual = Z;
yResidual = zY;
previousEnergy = real(yResidual' * yResidual);
identity = eye(nFeatures);

for k = 1:OMax
    correlation = ZResidual' * yResidual;
    if mu == 0
        w = correlation;
    else
        visibleScale = real(trace(bandGram)) / nFeatures;
        penalty = (fullGram + fullGram') / (2 * visibleScale);
        w = (identity + mu * penalty) \ correlation;
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

BByO = zeros(nFeatures, OMax);
for O = 1:OMax
    R = W(:, 1:O) / (P(:, 1:O)' * W(:, 1:O));
    T = Z * R;
    C = (T' * T + lambda * eye(O)) \ (T' * zY);
    BByO(:, O) = (R * C) ./ scaleA' * scaleY;
end
end

function O = select_adaptive_order(deltaList, kappaList, OMax, T_delta, kappaGMax)
% Follow the stopping rule in Adaptive latent variable selection/main.m.
firstRejected = find( ...
    (T_delta > 0 & deltaList(2:OMax) <= T_delta) | ...
    kappaList(2:OMax) >= kappaGMax, 1);
if isempty(firstRejected)
    O = OMax;
else
    O = firstRejected;
end
end
