function [y_LD_PLS, gain_rms, O] = get_adaptive_LD_PLS( ...
    K, M, X, x, h, y_bf, lambda, T_delta, kappa_G_max)
%GET_ADAPTIVE_LD_PLS 全带宽能量惩罚方向，自适应选择潜变量数。
% 惩罚强度固定为 mu=1。T_delta=0 关闭边际增益条件；任一生效的
% 停止条件触发时拒绝当前候选，并始终保留第一个潜变量。

H_x = geth_even_BLMP(x, K, M, h);
[hFull, O] = get_energy_regularized_LD_PLS( ...
    X, H_x, y_bf, lambda, 1, T_delta, kappa_G_max, size(X, 2) - 1);
y_LD_PLS = X * hFull;
gain_rms = rms(y_LD_PLS) / rms(x);

% 原无惩罚实现留作对照；正式入口不执行。实验需要无惩罚时可在共用函数中设 mu=0。
%{
xBandLinear = H_x(:, 1);
wLinear = xBandLinear \ y_bf;
bandYResidual = y_bf - xBandLinear * wLinear;

pBand = xBandLinear \ H_x(:, 2:end);
bandNonlinear = H_x(:, 2:end) - xBandLinear * pBand;
scaleX = std(bandNonlinear);
scaleY = std(bandYResidual);
Z = (bandNonlinear - mean(bandNonlinear)) ./ scaleX;
zY = (bandYResidual - mean(bandYResidual)) ./ scaleY;

maxO = size(Z, 2);
W = zeros(maxO);
P = zeros(maxO);
scores = zeros(size(Z, 1), maxO);
ZResidual = Z;
yResidual = zY;
previousEnergy = real(yResidual' * yResidual);
O = 0;

for k = 1:maxO
    w = ZResidual' * yResidual;
    w = w / norm(w);
    t = ZResidual * w;
    p = ZResidual' * t / (t' * t);
    q = yResidual' * t / (t' * t);
    newYResidual = yResidual - t * q;
    newEnergy = real(newYResidual' * newYResidual);
    Delta_k = (previousEnergy - newEnergy) / (previousEnergy + eps);

    candidateScores = [scores(:, 1:k-1), t];
    kappa_G_k = cond(candidateScores' * candidateScores);
    if k > 1 && (kappa_G_k >= kappa_G_max || ...
            (T_delta > 0 && Delta_k <= T_delta))
        break;
    end

    W(:, k) = w;
    P(:, k) = p;
    scores(:, k) = t;
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
B = (R * C) ./ scaleX' * scaleY;

pFull = X(:, 1) \ X(:, 2:end);
fullNonlinear = X(:, 2:end) - X(:, 1) * pFull;
y_LD_PLS = fullNonlinear * B + X(:, 1) * wLinear;
gain_rms = rms(y_LD_PLS) / rms(x);
%}
end
