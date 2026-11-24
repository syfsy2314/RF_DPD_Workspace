function fig = plot_psd_comparison(results)
%PLOT_PSD_COMPARISON 在可见 MATLAB figure 中绘制 PA 输出恢复 PSD。

    f = results.psdFrequencyMHz;
    psdDb = results.psdDb;
    narrowMHz = min(results.h2EdgesMHz);

    fig = figure('Name', 'PA 输出恢复 PSD', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1050, 600]);
    plot(f, psdDb(:, 1), 'k-', 'LineWidth', 1.2); hold on;
    plot(f, psdDb(:, 2), 'b-', 'LineWidth', 1.1);
    plot(f, psdDb(:, 3), 'r-', 'LineWidth', 1.1);
    title(sprintf('H2 %d MHz：PA 输出恢复', narrowMHz));
    legend('真实 y_0', '仅相位', '相位+幅度', 'Location', 'best');
    xlabel('频率 (MHz)');
    ylabel('PSD (dB/MHz)');
    xlim([-500, 500]);
    grid on;
end
