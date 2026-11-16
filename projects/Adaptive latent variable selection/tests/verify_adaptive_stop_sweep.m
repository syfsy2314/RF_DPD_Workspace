function study = verify_adaptive_stop_sweep()
% Calibrate shared stopping thresholds on validation data, then test once.

data = verify_pls_direction_regularization();
deltaThresholds = [0 logspace(-7, -2, 21)];
kappaThresholds = [1e2 logspace(2.5, 8, 23) inf];
nDelta = numel(deltaThresholds);
nKappa = numel(kappaThresholds);
nBandwidths = numel(data.bandwidths);
nMu = numel(data.muList);
maxO = data.OList(end);
originalValidation = zeros(nDelta, nKappa);
energyValidation = zeros(nDelta, nKappa);
originalValidationByBandwidth = zeros(nBandwidths, nDelta, nKappa);
energyValidationByBandwidth = zeros(nBandwidths, nDelta, nKappa);

for id = 1:nDelta
    for ik = 1:nKappa
        originalByBandwidth = zeros(nBandwidths, 1);
        energyByBandwidth = zeros(nBandwidths, 1);
        for ib = 1:nBandwidths
            validationByMu = zeros(1, nMu);
            for im = 1:nMu
                O = select_order( ...
                    data.deltaByStep(ib, :, im), data.kappaByStep(ib, :, im), ...
                    maxO, deltaThresholds(id), kappaThresholds(ik));
                validationByMu(im) = data.validationNmse(ib, O, im);
            end
            originalByBandwidth(ib) = validationByMu(1);
            energyByBandwidth(ib) = min(validationByMu);
        end
        originalValidation(id, ik) = mean(originalByBandwidth);
        energyValidation(id, ik) = mean(energyByBandwidth);
        originalValidationByBandwidth(:, id, ik) = originalByBandwidth;
        energyValidationByBandwidth(:, id, ik) = energyByBandwidth;
    end
end

[~, originalIndex] = min(originalValidation, [], 'all', 'linear');
[originalDeltaIndex, originalKappaIndex] = ind2sub([nDelta, nKappa], originalIndex);
[~, energyIndex] = min(energyValidation, [], 'all', 'linear');
[energyDeltaIndex, energyKappaIndex] = ind2sub([nDelta, nKappa], energyIndex);
originalThresholds = [deltaThresholds(originalDeltaIndex), kappaThresholds(originalKappaIndex)];
energyThresholds = [deltaThresholds(energyDeltaIndex), kappaThresholds(energyKappaIndex)];
[originalDeltaOnlyValidation, originalDeltaOnlyIndex] = min(originalValidation(:, end));
[energyDeltaOnlyValidation, energyDeltaOnlyIndex] = min(energyValidation(:, end));
leaveOutOriginalThresholds = zeros(nBandwidths, 2);
leaveOutEnergyThresholds = zeros(nBandwidths, 2);
leaveOutOriginalO = zeros(nBandwidths, 1);
leaveOutEnergyO = zeros(nBandwidths, 1);
leaveOutEnergyMu = zeros(nBandwidths, 1);
leaveOutOriginalTest = zeros(nBandwidths, 1);
leaveOutEnergyTest = zeros(nBandwidths, 1);
for heldOut = 1:nBandwidths
    tuningBandwidths = setdiff(1:nBandwidths, heldOut);
    originalScore = squeeze(mean(originalValidationByBandwidth(tuningBandwidths, :, :), 1));
    energyScore = squeeze(mean(energyValidationByBandwidth(tuningBandwidths, :, :), 1));
    [~, originalIndexHeldOut] = min(originalScore, [], 'all', 'linear');
    [~, energyIndexHeldOut] = min(energyScore, [], 'all', 'linear');
    [idOriginal, ikOriginal] = ind2sub([nDelta nKappa], originalIndexHeldOut);
    [idEnergy, ikEnergy] = ind2sub([nDelta nKappa], energyIndexHeldOut);
    leaveOutOriginalThresholds(heldOut, :) = ...
        [deltaThresholds(idOriginal) kappaThresholds(ikOriginal)];
    leaveOutEnergyThresholds(heldOut, :) = ...
        [deltaThresholds(idEnergy) kappaThresholds(ikEnergy)];
    leaveOutOriginalO(heldOut) = select_order( ...
        data.deltaByStep(heldOut, :, 1), data.kappaByStep(heldOut, :, 1), ...
        maxO, leaveOutOriginalThresholds(heldOut, 1), leaveOutOriginalThresholds(heldOut, 2));
    leaveOutOriginalTest(heldOut) = data.testNmse(heldOut, leaveOutOriginalO(heldOut), 1);
    validationByMu = zeros(1, nMu);
    ordersByMu = zeros(1, nMu);
    for im = 1:nMu
        ordersByMu(im) = select_order( ...
            data.deltaByStep(heldOut, :, im), data.kappaByStep(heldOut, :, im), ...
            maxO, leaveOutEnergyThresholds(heldOut, 1), leaveOutEnergyThresholds(heldOut, 2));
        validationByMu(im) = data.validationNmse(heldOut, ordersByMu(im), im);
    end
    [~, selectedMuIndex] = min(validationByMu);
    leaveOutEnergyO(heldOut) = ordersByMu(selectedMuIndex);
    leaveOutEnergyMu(heldOut) = data.muList(selectedMuIndex);
    leaveOutEnergyTest(heldOut) = data.testNmse( ...
        heldOut, leaveOutEnergyO(heldOut), selectedMuIndex);
end

originalO = zeros(nBandwidths, 1);
originalTest = zeros(nBandwidths, 1);
energyO = zeros(nBandwidths, 1);
energyMu = zeros(nBandwidths, 1);
energyTest = zeros(nBandwidths, 1);
for ib = 1:nBandwidths
    originalO(ib) = select_order( ...
        data.deltaByStep(ib, :, 1), data.kappaByStep(ib, :, 1), ...
        maxO, originalThresholds(1), originalThresholds(2));
    originalTest(ib) = data.testNmse(ib, originalO(ib), 1);

    validationByMu = zeros(1, nMu);
    ordersByMu = zeros(1, nMu);
    for im = 1:nMu
        ordersByMu(im) = select_order( ...
            data.deltaByStep(ib, :, im), data.kappaByStep(ib, :, im), ...
            maxO, energyThresholds(1), energyThresholds(2));
        validationByMu(im) = data.validationNmse(ib, ordersByMu(im), im);
    end
    [~, selectedMuIndex] = min(validationByMu);
    energyO(ib) = ordersByMu(selectedMuIndex);
    energyMu(ib) = data.muList(selectedMuIndex);
    energyTest(ib) = data.testNmse(ib, energyO(ib), selectedMuIndex);
end
coreO = verify_core_orders(data.bandwidths, originalThresholds);
assert(isequal(coreO, originalO), ...
    'The adaptive project function disagrees with the validation trajectory.');

fprintf('Shared original thresholds: Delta=%g, kappa=%g; mean validation=%g dB\n', ...
    originalThresholds, originalValidation(originalIndex));
fprintf('Shared energy thresholds:   Delta=%g, kappa=%g; mean validation=%g dB\n', ...
    energyThresholds, energyValidation(energyIndex));
fprintf('Delta-only validation: original %g dB (threshold %g), energy %g dB (threshold %g)\n', ...
    originalDeltaOnlyValidation, deltaThresholds(originalDeltaOnlyIndex), ...
    energyDeltaOnlyValidation, deltaThresholds(energyDeltaOnlyIndex));
fprintf('BW   original adaptive O/test   original scan O/test   energy adaptive O/mu/test   energy scan O/mu/test\n');
for ib = 1:nBandwidths
    fprintf('%3d      %2d / %8.2f           %2d / %8.2f         %2d / %5g / %8.2f       %2d / %5g / %8.2f\n', ...
        data.bandwidths(ib), originalO(ib), originalTest(ib), ...
        data.originalO(ib), data.originalTestNmse(ib), ...
        energyO(ib), energyMu(ib), energyTest(ib), ...
        data.selectedO(ib), data.selectedMu(ib), data.selectedTestNmse(ib));
end
fprintf('Leave-one-bandwidth-out: BW   original Delta/kappa/O/test   energy Delta/kappa/O/mu/test\n');
for ib = 1:nBandwidths
    fprintf('%3d   %6g / %8g / %2d / %8.2f      %6g / %8g / %2d / %5g / %8.2f\n', ...
        data.bandwidths(ib), leaveOutOriginalThresholds(ib, :), ...
        leaveOutOriginalO(ib), leaveOutOriginalTest(ib), ...
        leaveOutEnergyThresholds(ib, :), leaveOutEnergyO(ib), ...
        leaveOutEnergyMu(ib), leaveOutEnergyTest(ib));
end

study.deltaThresholds = deltaThresholds;
study.kappaThresholds = kappaThresholds;
study.originalValidation = originalValidation;
study.energyValidation = energyValidation;
study.originalValidationByBandwidth = originalValidationByBandwidth;
study.energyValidationByBandwidth = energyValidationByBandwidth;
study.originalThresholds = originalThresholds;
study.energyThresholds = energyThresholds;
study.originalDeltaOnlyValidation = originalDeltaOnlyValidation;
study.energyDeltaOnlyValidation = energyDeltaOnlyValidation;
study.originalO = originalO;
study.coreO = coreO;
study.originalTest = originalTest;
study.energyO = energyO;
study.energyMu = energyMu;
study.energyTest = energyTest;
study.leaveOutOriginalThresholds = leaveOutOriginalThresholds;
study.leaveOutEnergyThresholds = leaveOutEnergyThresholds;
study.leaveOutOriginalO = leaveOutOriginalO;
study.leaveOutEnergyO = leaveOutEnergyO;
study.leaveOutEnergyMu = leaveOutEnergyMu;
study.leaveOutOriginalTest = leaveOutOriginalTest;
study.leaveOutEnergyTest = leaveOutEnergyTest;
study.data = data;
end

function O = select_order(deltaList, kappaList, OMax, T_delta, kappaGMax)
firstRejected = find( ...
    (T_delta > 0 & deltaList(2:OMax) <= T_delta) | ...
    kappaList(2:OMax) >= kappaGMax, 1);
if isempty(firstRejected)
    O = OMax;
else
    O = firstRejected;
end
end

function coreO = verify_core_orders(bandwidths, thresholds)
projectDir = fileparts(fileparts(mfilename('fullpath')));
dataX = load(fullfile(projectDir, '..', 'LD-PLS', 'data', 'x.mat'), 'x');
dataY = load(fullfile(projectDir, '..', 'LD-PLS', 'data', 'y0.mat'), 'y0');
x = dataX.x(:);
y0 = dataY.y0(:);
K = 7;
M = 5;
X = geth_even_1d(x, K, M);
train = 1:floor(0.6 * numel(x));
coreO = zeros(numel(bandwidths), 1);
for ib = 1:numel(bandwidths)
    BW = bandwidths(ib);
    h1 = design_filter([BW + 5 BW + 10], 0.5, 1e-4);
    h2 = design_filter([BW BW + 3], 0.1, 1e-6);
    y_f = filter(h1, 1, y0);
    [y_f, ~, ~] = align_circular_signal(y_f, x);
    y_bf = filter(h2, 1, y_f);
    [~, ~, coreO(ib)] = get_adaptive_LD_PLS( ...
        K, M, X(train, :), x(train), h2, y_bf(train), 10, ...
        thresholds(1), thresholds(2));
end
end

function h = design_filter(edges, rippleDb, stopRipple)
fs = 1000;
passRipple = (10^(rippleDb / 20) - 1) / (10^(rippleDb / 20) + 1);
[order, frequency, amplitude, weights] = firpmord( ...
    edges / (fs / 2), [1 0], [passRipple stopRipple]);
h = firpm(order, frequency, amplitude, weights, {20});
end
