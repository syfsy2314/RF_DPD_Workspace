function [hFull, O] = get_energy_regularized_LD_PLS( ...
    X, H_x, y_bf, lambda, mu, T_delta, kappa_G_max, maxO)
%GET_ENERGY_REGULARIZED_LD_PLS 拟合全带宽能量惩罚方向的自适应 LD-PLS。
% X/H_x 是训练段的全带宽/经 H2 滤波的基函数，y_bf 是带限观测。
% lambda 为最终回归的岭系数，mu 为方向惩罚强度；mu=0 退化为原方向。
% 返回全带宽系数 hFull（X*hFull 为预测输出）及保留的潜变量数 O。
% 提取到 maxO 或候选潜变量触发停止准则时结束，触发的候选不计入 O。

% 先分离线性列，避免高能量线性项主导非线性潜变量的提取。
xBandLinear = H_x(:, 1);
wLinear = xBandLinear \ y_bf;
pBand = xBandLinear \ H_x(:, 2:end);
pFull = X(:, 1) \ X(:, 2:end);
bandNonlinear = H_x(:, 2:end) - xBandLinear * pBand;
fullNonlinear = X(:, 2:end) - X(:, 1) * pFull;
yBandResidual = y_bf - xBandLinear * wLinear;

scaleA = std(bandNonlinear);
scaleY = std(yBandResidual);
Z = (bandNonlinear - mean(bandNonlinear)) ./ scaleA;
zY = (yBandResidual - mean(yBandResidual)) ./ scaleY;
U = (fullNonlinear - mean(fullNonlinear)) ./ scaleA;
% 两个 Gram 矩阵分别描述全带宽能量和带限可观测能量。
fullGram = U' * U;
bandGram = Z' * Z;

nFeatures = size(Z, 2);
W = zeros(nFeatures, maxO);
P = zeros(nFeatures, maxO);
scores = zeros(size(Z, 1), maxO);
ZResidual = Z;
yResidual = zY;
previousEnergy = real(zY' * zY);
identity = eye(nFeatures);
O = 0;

for k = 1:maxO
    correlation = ZResidual' * yResidual;
    if mu == 0
        w = correlation;
    else
        % 用当前带限 Gram 的平均能量归一化惩罚。
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
    Delta_k = (previousEnergy - newEnergy) / (previousEnergy + eps);
    candidateScores = [scores(:, 1:k-1), t];
    kappa_G_k = cond(candidateScores' * candidateScores);
    % OR 准则：增益过小或条件数过大即拒绝当前潜变量；至少保留第一个。
    if k > 1 && (kappa_G_k >= kappa_G_max || ...
            (T_delta > 0 && Delta_k <= T_delta))
        break;
    end

    W(:, k) = w;
    P(:, k) = p;
    scores(:, k) = t;
    % PLS 消去当前方向后，同步更新两种能量度量。
    transformation = identity - w * p';
    fullGram = transformation' * fullGram * transformation;
    bandGram = transformation' * bandGram * transformation;
    ZResidual = ZResidual - t * p';
    yResidual = newYResidual;
    previousEnergy = newEnergy;
    O = k;
end

W = W(:, 1:O);
P = P(:, 1:O);
R = W / (P' * W);
T = Z * R;
C = (T' * T + lambda * eye(O)) \ (T' * zY);
B = (R * C) ./ scaleA' * scaleY;
% 将非线性残差系数还原到原始 X 的列空间，供测试段直接预测。
hFull = [wLinear - pFull * B; B];
end
