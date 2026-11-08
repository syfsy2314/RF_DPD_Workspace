function [yResults, gainResults] = run_amplitude_alternating_iterations( ...
        iterationCount, yInitial, gainInitial, yFeedback, h2, ...
        x, K, M, X, S, lamda, operatorR, operatorEta)
%RUN_AMPLITUDE_ALTERNATING_ITERATIONS 幅度算子交替迭代。
%   第 1 列为初始结果，第 i+1 列为第 i 次迭代结果。

    yResults = zeros(length(yInitial), iterationCount + 1, 'like', yInitial);
    gainResults = zeros(iterationCount + 1, 1);
    yResults(:, 1) = yInitial;
    gainResults(1) = gainInitial;

    for i = 1:iterationCount
        hModel = estimate_amplitude_operator( ...
            yResults(:, i), yFeedback, h2, operatorR, operatorEta);

        linearReference = filter(hModel, 1, x);
        yFeedbackAligned = align_circular_signal( ...
            yFeedback, linearReference);

        [yResults(:, i + 1), gainResults(i + 1)] = get_OPLS( ...
            K, M, X, x, hModel, yFeedbackAligned, S, lamda);
    end
end
