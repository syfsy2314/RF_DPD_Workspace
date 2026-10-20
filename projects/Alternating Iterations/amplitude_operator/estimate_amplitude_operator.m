function [hModel, hA, a, fitNmseDb] = ...
        estimate_amplitude_operator(yCoarse, yFeedback, h2, R, eta)
%ESTIMATE_AMPLITUDE_OPERATOR 由粗略 PA 输出估计实对称幅度算子。
%   yCoarse   ：第 0 轮粗略恢复的 PA 全带输出
%   yFeedback ：实际带限反馈信号
%   h2        ：已知数字覆盖滤波器
%   R         ：对称算子单边阶数，完整 FIR 长度为 2*R+1
%   eta       ：无量纲岭正则化系数

    yCoarse = yCoarse(:);
    yFeedback = yFeedback(:);
    h2 = h2(:);
    if numel(yCoarse) ~= numel(yFeedback)
        error('estimate_amplitude_operator:SizeMismatch', ...
            '粗略输出与反馈信号的样本数必须一致。');
    end

    %% 1. 构造实对称幅度算子的回归矩阵 Z
    % 先经过已知 H2，再使用正负对称移位表示零相位幅度响应。
    filteredOutput = filter(h2, 1, yCoarse);
    Z = zeros(numel(yCoarse), R + 1, 'like', yCoarse);
    Z(:, 1) = filteredOutput;
    for r = 1:R
        Z(:, r + 1) = circshift(filteredOutput, r) + ...
            circshift(filteredOutput, -r);
    end

    % 反馈信号只进行一次公共时延和常相位对齐，禁止逐列对齐 Z。
    [yTarget, ~, ~] = align_circular_signal(yFeedback, Z(:, 1));

    %% 2. 在有效样本区间进行一次实参数岭回归
    % 去除 FIR 起始暂态和 circshift 在记录两端产生的环绕样本。
    guardSamples = (numel(h2) - 1) + R;
    validIndices = (guardSamples + 1):(numel(yCoarse) - guardSamples);
    if isempty(validIndices)
        error('estimate_amplitude_operator:InsufficientSamples', ...
            '当前信号长度不足以覆盖所选算子的保护区间。');
    end

    ZValid = Z(validIndices, :);
    yValid = yTarget(validIndices);
    gram = real(ZValid' * ZValid);
    gram = (gram + gram') / 2;
    crossTerm = real(ZValid' * yValid);
    lambda = eta * trace(gram) / size(gram, 1);
    aRaw = real((gram + lambda * eye(size(gram), 'like', gram)) ...
        \ crossTerm);

    residual = yValid - ZValid * aRaw;
    fitNmseDb = 10 * log10(max( ...
        norm(residual)^2 / norm(yValid)^2, realmin));

    %% 3. 消除尺度歧义并生成模型侧复合滤波器
    % A(0)=a0+2*sum(ar)=1；hA 为长度 2*R+1 的因果偶对称 FIR。
    centerGain = aRaw(1) + 2 * sum(aRaw(2:end));
    if abs(centerGain) <= 1e-8 * max(1, norm(aRaw, 1))
        error('estimate_amplitude_operator:InvalidCenterGain', ...
            '估计算子的中心频率增益过小，无法完成归一化。');
    end
    a = aRaw / centerGain;
    hA = [flipud(a(2:end)); a];
    hModel = conv(hA, h2);
end
