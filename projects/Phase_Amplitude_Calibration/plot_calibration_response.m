function fig = plot_calibration_response(results)
%PLOT_CALIBRATION_RESPONSE 绘制 H1 全通带真实响应、校准观测与模型响应。

    fH1 = results.h1FrequencyMHz;
    fObserved = results.frequencyMHz;
    fig = figure('Name', 'H1 幅度与相位校准', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1050, 750]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');

    nexttile;
    plot(fH1, results.h1TrueMagnitudeDb, 'k-', ...
        'LineWidth', 1.2); hold on;
    plot(fObserved, results.magnitudeDb, 'b-', ...
        'LineWidth', 1.2);
    plot(fH1, results.h1ModelMagnitudeDb, 'r--', ...
        'LineWidth', 1.1);
    xline(min(fObserved), ':', 'Color', [0.5 0.5 0.5]);
    xline(max(fObserved), ':', 'Color', [0.5 0.5 0.5]);
    xlim([-85, 85]);
    ylabel('归一化幅度 (dB)');
    title('H1 通带幅度');
    legend('真实 H1', '开关 2 观测', '模型侧 FIR', ...
        'Location', 'best');
    grid on;

    nexttile;
    plot(fH1, results.h1TruePhaseResidualRad * 180 / pi, ...
        'k-', 'LineWidth', 1.2); hold on;
    plot(fObserved, results.phaseResidualRad * 180 / pi, ...
        'b-', 'LineWidth', 1.2);
    plot(fH1, results.h1ModelPhaseResidualRad * 180 / pi, ...
        'r--', 'LineWidth', 1.1);
    xline(min(fObserved), ':', 'Color', [0.5 0.5 0.5]);
    xline(max(fObserved), ':', 'Color', [0.5 0.5 0.5]);
    xlim([-85, 85]);
    xlabel('频率 (MHz)');
    ylabel('残余相位 (度)');
    title('H1 通带相位');
    grid on;
end
