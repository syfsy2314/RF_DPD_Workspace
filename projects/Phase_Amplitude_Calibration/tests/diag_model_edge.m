function metrics = diag_model_edge()
%DIAG_MODEL_EDGE 比较局部校准 FIR 与真实 H1 在 H2 通带边缘的响应。

    projectDir = fileparts(fileparts(mfilename('fullpath')));
    d = load(fullfile(projectDir, 'results', 'results.mat'));
    r = d.results;
    Fs = 1000;
    nfft = 8192;
    tapCount = 4001;
    taperMHz = 5;
    f = r.frequencyMHz;
    design = [ones(numel(f), 1), 2 * pi * f / Fs];
    hCombined = create_calibration_fir(f, r.phaseResidualRad, ...
        r.magnitudeDb, Fs, nfft, tapCount, taperMHz);
    H = freqz(hCombined, 1, f, Fs);
    modelPhase = unwrap(angle(H));
    modelPhase = modelPhase - design * (design \ modelPhase);
    modelMagnitudeDb = 20 * log10(abs(H));
    [~, centerIndex] = min(abs(f));
    modelMagnitudeDb = modelMagnitudeDb - modelMagnitudeDb(centerIndex);
    edgeMHz = min(r.h2EdgesMHz);
    interiorMHz = edgeMHz - 5;
    masks = {abs(f) <= interiorMHz, ...
        abs(f) > interiorMHz & abs(f) <= edgeMHz};
    names = {'interior', 'edge'};
    for i = 1:2
        mask = masks{i};
        metrics.([names{i}, 'PhaseErrorRad']) = rms( ...
            modelPhase(mask) - r.truePhaseResidualRad(mask));
        metrics.([names{i}, 'MagnitudeErrorDb']) = rms( ...
            modelMagnitudeDb(mask) - r.trueMagnitudeDb(mask));
    end
    disp(metrics);
end
