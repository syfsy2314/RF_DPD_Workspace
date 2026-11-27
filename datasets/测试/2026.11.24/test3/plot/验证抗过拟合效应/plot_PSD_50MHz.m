%% Load the paired IQ and 50 MHz recovery results
% 原始输入 x 来自 datasets/data_38400_v0.6_test1.dat，全带 PA 参考 y0 来自 test2/y0.mat。
% MG、MR、LD-PLS、SBL 的 50 MHz 信号基于 test2 配对 IQ 离线恢复。
% 恢复 MAT 文件归档在 test3 各方法的 recovery 目录里，变量 y 是离线恢复结果。
scriptDir = fileparts(mfilename('fullpath'));
testDir = fileparts(fileparts(scriptDir));
sessionDir = fileparts(testDir);
datasetsDir = fileparts(fileparts(sessionDir));
run(fullfile(fileparts(datasetsDir), 'startup.m'));

reference = load(fullfile(sessionDir, 'test2', 'y0.mat'), 'y0');
y0 = reference.y0(:);
[~, ~, x] = import_iq_binary(fullfile(datasetsDir, 'data_38400_v0.6_test1.dat'), numel(y0));
x = x / max(abs(x));

recovered = load(fullfile(testDir, 'MG', 'recovery', ...
    'y_test2_bw50_mg_s3.mat'), 'y');
y_MG = recovered.y(:);
recovered = load(fullfile(testDir, 'MR', 'recovery', ...
    'y_test2_bw50_mr_a1_0p0001_a2_0.mat'), 'y');
y_MR = recovered.y(:);
recovered = load(fullfile(testDir, 'LD-PLS', 'recovery', ...
    'y_test2_bw50_ldpls_mu0p01_O8.mat'), 'y');
y_OPLS = recovered.y(:);
recovered = load(fullfile(testDir, 'SBL', 'recovery', ...
    'y_test2_bw50_sbl_iter3.mat'), 'y');
y_SBL = recovered.y(:);

fs = 1000;
nfft = 1024;
window=hann(nfft);
N = length(x);
gain_rmsy0 = rms(y0)/rms(x);
gain_rmsMG = rms(y_MG)/rms(x);
gain_rmsMR = rms(y_MR)/rms(x);
gain_rmsOPLS = rms(y_OPLS)/rms(x);
gain_rmsSBL = rms(y_SBL)/rms(x);
%%
figure(1);
set(gcf, 'Color', 'w');
[p3, f3] = pwelch(y_MG/gain_rmsMG, window, nfft/2, nfft ,fs);
plot(f3 - 1/2 * fs, 10 * log10(fftshift(p3*100)), 'c-')
hold on;
[p6, f6] = pwelch(y_OPLS/gain_rmsOPLS, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'm-')
hold on;
[p6, f6] = pwelch(y_SBL/gain_rmsSBL, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'g-')
hold on;
[p6, f6] = pwelch(y_MR/gain_rmsMR, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'b-')
hold on;
[p2, f2] = pwelch(y0/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f2 - 1/2 * fs, 10 * log10(fftshift(p2*100)), 'r-')
hold on;
xlabel('Frequency Offset(MHz)');
ylabel('PSD(dBm/Hz)');
legend('MG','LD-PLS','SBL','MR','FS Ref');
grid on;
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.15);
axis([-300 300 -100 20]);
set(gca, 'Color', 'w');
%%
figure(2);
set(gcf, 'Color', 'w');

[p3, f3] = pwelch(y_MG/gain_rmsMG-y0/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f3 - 1/2 * fs, 10 * log10(fftshift(p3*100)), 'c-')
hold on;
[p6, f6] = pwelch(y_OPLS/gain_rmsOPLS-y0/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'm-')
hold on;
[p6, f6] = pwelch(y_SBL/gain_rmsSBL-y0/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'g-')
hold on;
[p6, f6] = pwelch(y_MR/gain_rmsMR-y0/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'b-')
hold on;
[p2, f2] = pwelch(y0/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f2 - 1/2 * fs, 10 * log10(fftshift(p2*100)), 'r-')
hold on;

xlabel('Frequency Offset(MHz)');
ylabel('PSD(dBm/Hz)');
legend('MG','LD-PLS','SBL','MR','FS Ref');
grid on;
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.15);
axis([-300 300 -80 20]);
set(gca, 'Color', 'w');