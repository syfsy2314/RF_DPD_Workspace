function [h, fitNmseDb] = estimate_full_fir_operator(inputSignal, outputSignal, L, eta, validIndices)
%ESTIMATE_FULL_FIR_OPERATOR 由已知输入输出拟合完整复数 FIR 算子。
%   使用线性卷积和岭正则，避开显式构造样本数乘抽头数的回归矩阵。

    inputSignal = inputSignal(:);
    outputSignal = outputSignal(:);
    N = numel(inputSignal);
    nfft = 2^nextpow2(N + L - 1);
    inputSpectrum = fft(inputSignal, nfft);
    lambda = eta * norm(inputSignal)^2;
    target = [outputSignal(validIndices); zeros(L, 1)];
    operator = @(v, direction) applyOperator( ...
        v, direction, inputSpectrum, N, L, nfft, validIndices, lambda);
    [h, flag] = lsqr(operator, target, 1e-8, 300);
    if flag ~= 0 && flag ~= 1
        warning('estimate_full_fir_operator:LsqrFlag', ...
            '完整 FIR 拟合未充分收敛，lsqr 标志为 %d。', flag);
    end
    fitted = filter(h, 1, inputSignal);
    fitNmseDb = 10 * log10(norm(fitted(validIndices) - ...
        outputSignal(validIndices))^2 / norm(outputSignal(validIndices))^2);
end

function value = applyOperator(v, direction, inputSpectrum, N, L, nfft, validIndices, lambda)
    if strcmp(direction, 'notransp')
        hPadded = [v; zeros(nfft - L, 1)];
        filtered = ifft(inputSpectrum .* fft(hPadded));
        value = [filtered(validIndices); sqrt(lambda) * v];
    else
        residual = zeros(N, 1, 'like', v);
        residual(validIndices) = v(1:numel(validIndices));
        correlation = ifft(conj(inputSpectrum) .* fft(residual, nfft));
        value = correlation(1:L) + sqrt(lambda) * v(numel(validIndices) + 1:end);
    end
end
