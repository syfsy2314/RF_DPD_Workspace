function [y_LD_PLS, gain_rms] = get_LD_PLS( ...
    K, M, X, x, h, y_bf, S, lambda)
%GET_LD_PLS 固定潜变量数的全带宽能量惩罚 LD-PLS。
% O = round(S * K * (M + 1))；方向惩罚强度固定为 mu=1。

H_x = geth_even_BLMP(x, K, M, h);
O = round(S * size(H_x, 2));
[hFull, ~] = get_energy_regularized_LD_PLS( ...
    X, H_x, y_bf, lambda, 1, 0, inf, O);
y_LD_PLS = X * hFull;
gain_rms = rms(y_LD_PLS) / rms(x);

% 原无惩罚实现留作对照；正式入口不执行。实验需要无惩罚时可在共用函数中设 mu=0。
%{
xBandLinear = H_x(:, 1);
wLinear = xBandLinear \ y_bf;
yResidual = y_bf - xBandLinear * wLinear;
pBand = xBandLinear \ H_x(:, 2:end);
bandNonlinear = H_x(:, 2:end) - xBandLinear * pBand;
B = pls_B(bandNonlinear, yResidual, O, lambda);
pFull = X(:, 1) \ X(:, 2:end);
fullNonlinear = X(:, 2:end) - X(:, 1) * pFull;
y_LD_PLS = fullNonlinear * B + X(:, 1) * wLinear;
%}
end
