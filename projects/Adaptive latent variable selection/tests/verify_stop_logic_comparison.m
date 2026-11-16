function comparison = verify_stop_logic_comparison()
% Compare joint and either-condition stopping on identical PLS trajectories.
% Tune shared thresholds on validation data and evaluate the test data once.

data = verify_pls_direction_regularization();
deltaThresholds = [0 logspace(-7, 0, 29)];
kappaThresholds = [1 logspace(0.5, 8, 31) inf];
rules = {'both', 'either'};
fixedThresholds = [9e-4 1e3];
nBandwidths = numel(data.bandwidths);
nMu = numel(data.muList);
nDelta = numel(deltaThresholds);
nKappa = numel(kappaThresholds);
maxO = data.OList(end);

comparison.bandwidths = data.bandwidths;
comparison.deltaThresholds = deltaThresholds;
comparison.kappaThresholds = kappaThresholds;
comparison.fixedThresholds = fixedThresholds;

for ir = 1:numel(rules)
    rule = rules{ir};
    originalValidation = zeros(nDelta, nKappa);
    energyValidation = zeros(nDelta, nKappa);
    for id = 1:nDelta
        for ik = 1:nKappa
            originalByBandwidth = zeros(nBandwidths, 1);
            energyByBandwidth = zeros(nBandwidths, 1);
            for ib = 1:nBandwidths
                [originalByBandwidth(ib), energyByBandwidth(ib)] = ...
                    validation_scores(data, ib, maxO, ...
                    deltaThresholds(id), kappaThresholds(ik), rule, nMu);
            end
            originalValidation(id, ik) = mean(originalByBandwidth);
            energyValidation(id, ik) = mean(energyByBandwidth);
        end
    end

    [originalBestValidation, originalIndex] = min(originalValidation, [], 'all', 'linear');
    [energyBestValidation, energyIndex] = min(energyValidation, [], 'all', 'linear');
    [originalDeltaIndex, originalKappaIndex] = ind2sub( ...
        [nDelta nKappa], originalIndex);
    [energyDeltaIndex, energyKappaIndex] = ind2sub( ...
        [nDelta nKappa], energyIndex);
    originalThresholds = [deltaThresholds(originalDeltaIndex) ...
        kappaThresholds(originalKappaIndex)];
    energyThresholds = [deltaThresholds(energyDeltaIndex) ...
        kappaThresholds(energyKappaIndex)];

    ruleResult.fixed = evaluate_thresholds(data, fixedThresholds, ...
        fixedThresholds, rule, maxO, nMu);
    ruleResult.tuned = evaluate_thresholds(data, originalThresholds, ...
        energyThresholds, rule, maxO, nMu);
    ruleResult.tuned.originalThresholds = originalThresholds;
    ruleResult.tuned.energyThresholds = energyThresholds;
    ruleResult.tuned.originalMeanValidation = originalBestValidation;
    ruleResult.tuned.energyMeanValidation = energyBestValidation;
    ruleResult.originalValidation = originalValidation;
    ruleResult.energyValidation = energyValidation;
    comparison.(rule) = ruleResult;
end

fprintf('Fixed thresholds: Delta=%g, kappa=%g\n', fixedThresholds);
print_results(comparison, data, 'fixed');
fprintf('Tuned shared thresholds on validation data:\n');
for ir = 1:numel(rules)
    rule = rules{ir};
    tuned = comparison.(rule).tuned;
    fprintf('%6s original Delta=%g, kappa=%g, validation=%.2f dB; ', ...
        rule, tuned.originalThresholds, tuned.originalMeanValidation);
    fprintf('energy Delta=%g, kappa=%g, validation=%.2f dB\n', ...
        tuned.energyThresholds, tuned.energyMeanValidation);
end
print_results(comparison, data, 'tuned');
end

function [originalScore, energyScore] = validation_scores( ...
    data, ib, maxO, T_delta, kappaGMax, rule, nMu)
validationByMu = zeros(1, nMu);
for im = 1:nMu
    O = select_order(data.deltaByStep(ib, :, im), ...
        data.kappaByStep(ib, :, im), maxO, T_delta, kappaGMax, rule);
    validationByMu(im) = data.validationNmse(ib, O, im);
end
originalScore = validationByMu(1);
energyScore = min(validationByMu);
end

function result = evaluate_thresholds(data, originalThresholds, ...
    energyThresholds, rule, maxO, nMu)
nBandwidths = numel(data.bandwidths);
result.originalO = zeros(nBandwidths, 1);
result.originalValidation = zeros(nBandwidths, 1);
result.originalTest = zeros(nBandwidths, 1);
result.energyO = zeros(nBandwidths, 1);
result.energyMu = zeros(nBandwidths, 1);
result.energyValidation = zeros(nBandwidths, 1);
result.energyTest = zeros(nBandwidths, 1);
for ib = 1:nBandwidths
    O = select_order(data.deltaByStep(ib, :, 1), ...
        data.kappaByStep(ib, :, 1), maxO, ...
        originalThresholds(1), originalThresholds(2), rule);
    result.originalO(ib) = O;
    result.originalValidation(ib) = data.validationNmse(ib, O, 1);
    result.originalTest(ib) = data.testNmse(ib, O, 1);

    ordersByMu = zeros(1, nMu);
    validationByMu = zeros(1, nMu);
    for im = 1:nMu
        ordersByMu(im) = select_order(data.deltaByStep(ib, :, im), ...
            data.kappaByStep(ib, :, im), maxO, ...
            energyThresholds(1), energyThresholds(2), rule);
        validationByMu(im) = data.validationNmse(ib, ordersByMu(im), im);
    end
    [result.energyValidation(ib), selectedMuIndex] = min(validationByMu);
    result.energyO(ib) = ordersByMu(selectedMuIndex);
    result.energyMu(ib) = data.muList(selectedMuIndex);
    result.energyTest(ib) = data.testNmse( ...
        ib, result.energyO(ib), selectedMuIndex);
end
end

function O = select_order(deltaList, kappaList, maxO, T_delta, kappaGMax, rule)
smallGain = T_delta > 0 & deltaList(2:maxO) <= T_delta;
poorCondition = kappaList(2:maxO) >= kappaGMax;
if strcmp(rule, 'both')
    rejected = smallGain & poorCondition;
else
    rejected = smallGain | poorCondition;
end
firstRejected = find(rejected, 1);
if isempty(firstRejected)
    O = maxO;
else
    O = firstRejected;
end
end

function print_results(comparison, data, selection)
fprintf('%s: rule  BW   original O / test dB   energy O / mu / test dB\n', selection);
rules = {'both', 'either'};
for ir = 1:numel(rules)
    rule = rules{ir};
    result = comparison.(rule).(selection);
    for ib = 1:numel(data.bandwidths)
        fprintf('       %6s %3d       %2d / %8.2f       %2d / %5g / %8.2f\n', ...
            rule, data.bandwidths(ib), result.originalO(ib), ...
            result.originalTest(ib), result.energyO(ib), ...
            result.energyMu(ib), result.energyTest(ib));
    end
end
end
