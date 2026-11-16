function [y_LD_PLS, gain_rms] = get_LD_PLS( ...
    K, M, X, x, h, y_bf, S, lambda)
%GET_LD_PLS Fit LD-PLS with O = round(S * K * (M + 1)) latent variables.

H_x = geth_even_BLMP(x, K, M, h);
xBandLinear = H_x(:, 1);
wLinear = xBandLinear \ y_bf;
yResidual = y_bf - xBandLinear * wLinear;

pBand = xBandLinear \ H_x(:, 2:end);
bandNonlinear = H_x(:, 2:end) - xBandLinear * pBand;
O = round(S * size(H_x, 2));
B = pls_B(bandNonlinear, yResidual, O, lambda);

pFull = X(:, 1) \ X(:, 2:end);
fullNonlinear = X(:, 2:end) - X(:, 1) * pFull;
y_LD_PLS = fullNonlinear * B + X(:, 1) * wLinear;
gain_rms = rms(y_LD_PLS) / rms(x);
end
