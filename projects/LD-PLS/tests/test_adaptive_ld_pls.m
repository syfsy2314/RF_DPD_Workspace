function tests = test_adaptive_ld_pls
tests = functiontests(localfunctions);
end

function testZeroGainThresholdUsesConditionNumberOnly(testCase)
[K, M, X, x, y] = example_signals();
lambda = 10;
[yAdaptive, gainAdaptive, O] = get_adaptive_LD_PLS( ...
    K, M, X, x, 1, y, lambda, 0, inf);
[yFixed, gainFixed] = get_LD_PLS( ...
    K, M, X, x, 1, y, O / size(X, 2), lambda);

verifyEqual(testCase, O, size(X, 2) - 1);
verifySize(testCase, yAdaptive, size(x));
verifyLessThan(testCase, norm(yAdaptive - yFixed) / norm(yFixed), 1e-9);
verifyEqual(testCase, gainAdaptive, gainFixed, 'RelTol', 1e-9);
end

function testEitherThresholdRejectsSecondCandidate(testCase)
[K, M, X, x, y] = example_signals();
[~, ~, byCondition] = get_adaptive_LD_PLS( ...
    K, M, X, x, 1, y, 10, 0, 1);
[~, ~, byGain] = get_adaptive_LD_PLS( ...
    K, M, X, x, 1, y, 10, 1, inf);

verifyEqual(testCase, byCondition, 1);
verifyEqual(testCase, byGain, 1);
end

function [K, M, X, x, y] = example_signals()
rng(17);
K = 3;
M = 1;
x = randn(512, 1) + 1i * randn(512, 1);
X = geth_even_1d(x, K, M);
y = X * [1; 0.3; -0.1; 0.2; 0.05; -0.03] + ...
    0.01 * (randn(size(x)) + 1i * randn(size(x)));
end
