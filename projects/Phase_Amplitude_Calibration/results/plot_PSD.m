% 从已有结果打开可交互的 MATLAB PSD figure，无需重跑模型辨识。
resultDir = fileparts(mfilename('fullpath'));
data = load(fullfile(resultDir, 'results.mat'));
plot_psd_comparison(data.results);
