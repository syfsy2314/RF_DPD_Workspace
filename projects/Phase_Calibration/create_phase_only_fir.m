function [hPhase, magnitudeDeviationDb] = create_phase_only_fir( ...
        fReliable, phaseResidual, Fs, nfft, tapCount, taperMHz)
%CREATE_PHASE_ONLY_FIR 将可靠频段的残余相位转换为因果 FIR 算子。
%   未观测频段的相位修正为零，可靠频段边缘使用线性过渡。

    fLimit = min(-min(fReliable), max(fReliable));
    fGrid = (-nfft / 2:nfft / 2 - 1).' * Fs / nfft;
    phaseGrid = zeros(nfft, 1);
    inside = abs(fGrid) <= fLimit;
    phaseGrid(inside) = interp1(fReliable, phaseResidual, ...
        fGrid(inside), 'pchip');
    taper = min(1, max(0, (fLimit - abs(fGrid)) / taperMHz));
    phaseSpectrum = exp(1j * phaseGrid .* taper);
    centeredImpulse = fftshift(ifft(ifftshift(phaseSpectrum)));
    center = nfft / 2 + 1;
    halfTaps = (tapCount - 1) / 2;
    hPhase = centeredImpulse(center - halfTaps:center + halfTaps);

    response = freqz(hPhase, 1, fReliable, Fs);
    magnitudeDeviationDb = max(abs(20 * log10(abs(response))));
end
