function tests = test_adaptive_stop
tests = functiontests(localfunctions);
end

function testEitherThresholdStopsAfterFirstScore(testCase)
[K, M, X, x, y] = example_signals();
S = 0.5;
lambda = 1;

[~, ~, byCondition] = get_ALVS_LD_PLS( ...
    K, M, X, x, 1, y, S, lambda, 0, 0);
[~, ~, byGain] = get_ALVS_LD_PLS( ...
    K, M, X, x, 1, y, S, lambda, 1, inf);

verifyEqual(testCase, byCondition, 1);
verifyEqual(testCase, byGain, 1);
end

function testDisabledThresholdsReachMaximum(testCase)
[K, M, X, x, y] = example_signals();
S = 0.5;
[~, ~, selected] = get_ALVS_LD_PLS( ...
    K, M, X, x, 1, y, S, 1, 0, inf);

verifyEqual(testCase, selected, round(S * size(X, 2)));
end

function [K, M, X, x, y] = example_signals()
rng(7);
K = 3;
M = 1;
x = randn(256, 1) + 1i * randn(256, 1);
X = geth_even_1d(x, K, M);
y = X * [1; 0.2; 0.05; 0.3; -0.1; 0.02];
end
