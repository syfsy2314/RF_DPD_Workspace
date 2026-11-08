function [hModel, hA, a, fitNmseDb] = ...
        estimate_complex_operator(yCoarse, yFeedback, h2, R, eta)
%ESTIMATE_COMPLEX_OPERATOR 由粗略 PA 输出盲估完整算子（幅度 + 相位）。
%  与 estimate_amplitude_operator 的区别：移位基不强制对称，保留相位信息，
%  因此能同时补偿 IIR 模拟滤波器的非线性相位（群延迟变化）。
%  不需要 H1 的任何先验信息，完全从信号盲估。
%   yCoarse   ：粗略恢复的 PA 全带输出（简单覆盖结果）
%   yFeedback ：实际带限反馈信号
%   h2        ：已知数字覆盖滤波器
%   R         ：算子单边阶数，完整 FIR 长度为 2*R+1
%   eta       ：无量纲岭正则化系数

    yCoarse = yCoarse(:);
    yFeedback = yFeedback(:);
    h2 = h2(:);
    if numel(yCoarse) ~= numel(yFeedback)
        error('estimate_complex_operator:SizeMismatch', ...
            '粗略输出与反馈信号的样本数必须一致。');
    end

    %% 1. 构造非对称移位回归矩阵 Z（覆盖 -R..R，保留相位）
    filteredOutput = filter(h2, 1, yCoarse);
    Z = zeros(numel(yCoarse), 2 * R + 1, 'like', yCoarse);
    for r = -R:R
        Z(:, r + R + 1) = circshift(filteredOutput, r);
    end

    % 反馈信号对齐到中心列（r=0 即 filteredOutput 本身）
    [yTarget, ~, ~] = align_circular_signal(yFeedback, Z(:, R + 1));

    %% 2. 复数岭回归（不取实部，保留相位）
    guardSamples = (numel(h2) - 1) + R;
    validIndices = (guardSamples + 1):(numel(yCoarse) - guardSamples);
    if isempty(validIndices)
        error('estimate_complex_operator:InsufficientSamples', ...
            '当前信号长度不足以覆盖所选算子的保护区间。');
    end

    ZValid = Z(validIndices, :);
    yValid = yTarget(validIndices);
    gram = ZValid' * ZValid;
    gram = (gram + gram') / 2;              % 确保 Hermitian
    crossTerm = ZValid' * yValid;
    lambda = eta * trace(real(gram)) / size(gram, 1);
    aRaw = (gram + lambda * eye(size(gram), 'like', gram)) \ crossTerm;

    residual = yValid - ZValid * aRaw;
    fitNmseDb = 10 * log10(max( ...
        norm(residual)^2 / norm(yValid)^2, realmin));

    %% 3. 尺度归一化（中心抽头），生成模型侧复合滤波器
    centerGain = aRaw(R + 1);
    if abs(centerGain) <= 1e-8 * max(1, norm(aRaw, 1))
        error('estimate_complex_operator:InvalidCenterGain', ...
            '估计算子的中心增益过小，无法完成归一化。');
    end
    a = aRaw / centerGain;
    hA = a;                                   % 非对称，含相位
    hModel = conv(hA, h2);
end
