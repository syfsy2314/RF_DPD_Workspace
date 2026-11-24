function h = create_calibration_fir(fReliable, phaseResidual, magnitudeDb, ...
        Fs, nfft, tapCount, taperMHz)
%CREATE_CALIBRATION_FIR 将开关 2 测得的局部幅度和残余相位构造成模型侧 FIR。
% 未观测频段取单位响应，可靠频段边缘逐渐退回单位响应。

    fLimit = min(-min(fReliable), max(fReliable));
    fGrid = (-nfft / 2:nfft / 2 - 1).' * Fs / nfft;
    phaseGrid = zeros(nfft, 1);
    magnitudeGridDb = zeros(nfft, 1);
    inside = abs(fGrid) <= fLimit;
    phaseGrid(inside) = interp1(fReliable, phaseResidual, ...
        fGrid(inside), 'pchip');
    magnitudeGridDb(inside) = interp1(fReliable, magnitudeDb, ...
        fGrid(inside), 'pchip');
    taper = min(1, max(0, (fLimit - abs(fGrid)) / taperMHz));
    spectrum = 10.^(magnitudeGridDb .* taper / 20) .* ...
        exp(1j * phaseGrid .* taper);
    centeredImpulse = fftshift(ifft(ifftshift(spectrum)));
    center = nfft / 2 + 1;
    halfTaps = (tapCount - 1) / 2;
    h = centeredImpulse(center - halfTaps:center + halfTaps);
end
