% 从已有结果打开 H1 全通带幅相图，无需重跑模型辨识。
resultDir = fileparts(mfilename('fullpath'));
data = load(fullfile(resultDir, 'results.mat'));
plot_calibration_response(data.results);
