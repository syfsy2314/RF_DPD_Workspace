function [yResults, gainResults, latentCounts] = run_amplitude_alternating_iterations( ...
        iterationCount, yInitial, gainInitial, yFeedback, h2, ...
        x, K, M, X, lambda, T_delta, kappa_G_max, operatorR, operatorEta)
%RUN_AMPLITUDE_ALTERNATING_ITERATIONS 幅度算子交替迭代。
%   第 1 列为初始结果，第 i+1 列为第 i 次迭代结果。

    yResults = zeros(length(yInitial), iterationCount + 1, 'like', yInitial);
    gainResults = zeros(iterationCount + 1, 1);
    latentCounts = zeros(iterationCount + 1, 1);
    latentCounts(1) = NaN;     % 初始估计的阶数由调用方记录
    yResults(:, 1) = yInitial;
    gainResults(1) = gainInitial;

    for i = 1:iterationCount
        hModel = estimate_amplitude_operator( ...
            yResults(:, i), yFeedback, h2, operatorR, operatorEta);

        linearReference = filter(hModel, 1, x);
        yFeedbackAligned = align_circular_signal( ...
            yFeedback, linearReference);

        [yResults(:, i + 1), gainResults(i + 1), latentCounts(i + 1)] = ...
            get_adaptive_LD_PLS(K, M, X, x, hModel, yFeedbackAligned, ...
            lambda, T_delta, kappa_G_max);
    end
end
