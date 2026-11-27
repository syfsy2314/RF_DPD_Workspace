clear;
close all;

%% 导入本次 50 MHz 实测 DPD 回采
% 原始输入 x 来自 datasets/data_38400_v0.6_test1.dat。
% 未做 DPD 的 PA 输出和全带 MP 基准分别来自 test2/MP/y0.bin、test2/MP/y0_dpd.bin。
% 四种方法的实测 DPD 结果来自 test3 各方法的 capture/y_dpd_50.bin。
scriptDir = fileparts(mfilename('fullpath'));
testDir = fileparts(fileparts(scriptDir));
sessionDir = fileparts(testDir);
datasetsDir = fileparts(fileparts(sessionDir));
run(fullfile(fileparts(datasetsDir), 'startup.m'));

sampleCount = 38400;
[~, ~, x] = import_iq_binary(fullfile(datasetsDir, 'data_38400_v0.6_test1.dat'), sampleCount);
x = x / max(abs(x));

mpDir = fullfile(sessionDir, 'test2', 'MP');
[~, ~, y0] = import_iq_binary(fullfile(mpDir, 'y0.bin'), sampleCount);
[~, ~, y0_DPD] = import_iq_binary(fullfile(mpDir, 'y0_dpd.bin'), sampleCount);
[~, ~, y_MG_DPD] = import_iq_binary(fullfile(testDir, 'MG', 'capture', 'y_dpd_50.bin'), sampleCount);
[~, ~, y_OPLS_DPD] = import_iq_binary(fullfile(testDir, 'LD-PLS', 'capture', 'y_dpd_50.bin'), sampleCount);
[~, ~, y_SBL_DPD] = import_iq_binary(fullfile(testDir, 'SBL', 'capture', 'y_dpd_50.bin'), sampleCount);
[~, ~, y_MR_DPD] = import_iq_binary(fullfile(testDir, 'MR', 'capture', 'y_dpd_50.bin'), sampleCount);

y0 = align_circular_signal(y0, x);
y0_DPD = align_circular_signal(y0_DPD, x);
y_MG_DPD = align_circular_signal(y_MG_DPD, x);
y_OPLS_DPD = align_circular_signal(y_OPLS_DPD, x);
y_SBL_DPD = align_circular_signal(y_SBL_DPD, x);
y_MR_DPD = align_circular_signal(y_MR_DPD, x);

gain_rmspa = rms(y0) / rms(x);
gain_rmsy0 = rms(y0_DPD) / rms(x);
gain_rmsMG = rms(y_MG_DPD) / rms(x);
gain_rmsOPLS = rms(y_OPLS_DPD) / rms(x);
gain_rmsSBL = rms(y_SBL_DPD) / rms(x);
gain_rmsMR = rms(y_MR_DPD) / rms(x);

%%
fs = 1000;
nfft = 1024;
window=hann(nfft);
N = length(x);
%%
ACPR1=comm.ACPR('SampleRate',fs*1e6,...
    'MainChannelFrequency',0,...
    'MainMeasurementBandwidth',100e6,...
    'AdjacentMeasurementBandwidth',100e6,...
    'AdjacentChannelOffset',[-100e6,100e6],...
    'MainChannelPowerOutputPort',true,...
    'AdjacentChannelPowerOutputPort',true);
ACPR_Original = step(ACPR1,x);
ACPR_NO_DPD = step(ACPR1,y0);
ACPR_ref_DPD = step(ACPR1,y0_DPD);
ACPR_MG_DPD = step(ACPR1,y_MG_DPD);
ACPR_LD_PLS_DPD = step(ACPR1,y_OPLS_DPD);
ACPR_SBL_DPD = step(ACPR1,y_SBL_DPD);
ACPR_MR_DPD = step(ACPR1,y_MR_DPD);
%%
figure(1);
set(gcf, 'Color', 'w');
[p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
hold on;
[p2, f2] = pwelch(y0_DPD/gain_rmsy0, window, nfft/2, nfft ,fs);
plot(f2 - 1/2 * fs, 10 * log10(fftshift(p2*100)), 'b-')
hold on;
[p3, f3] = pwelch(y_MG_DPD/gain_rmsMG, window, nfft/2, nfft ,fs);
plot(f3 - 1/2 * fs, 10 * log10(fftshift(p3*100)), 'c-')
hold on;
[p6, f6] = pwelch(y_OPLS_DPD/gain_rmsOPLS, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'm-')
hold on;
[p6, f6] = pwelch(y_SBL_DPD/gain_rmsSBL, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'g-')
hold on;
[p6, f6] = pwelch(y_MR_DPD/gain_rmsMR, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'k-')
hold on;
xlabel('Frequency Offset(MHz)');
ylabel('PSD(dBm/Hz)');
legend('W/o DPD','FS-DPD Ref','MG-DPD','LD-PLS-DPD','SBL-DPD','MR-DPD');
grid on;
set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.15);
axis([-300 300 -80 20]);
set(gca, 'Color', 'w');
%%
EVM_NO_DPD = get_EVM(y0,x);
EVM_ref_DPD = get_EVM(y0_DPD,x);
EVM_MG_DPD = get_EVM(y_MG_DPD,x);
EVM_LD_PLS_DPD = get_EVM(y_OPLS_DPD,x);
EVM_SBL_DPD = get_EVM(y_SBL_DPD,x);
EVM_MR_DPD = get_EVM(y_MR_DPD,x);
%% nmse
NMSE_MG_DPD = get_NMSE(y_MG_DPD,x);
NMSE_LD_PLS_DPD = get_NMSE(y_OPLS_DPD,x);
NMSE_SBL_DPD = get_NMSE(y_SBL_DPD,x);
NMSE_MR_DPD = get_NMSE(y_MR_DPD,x);
NMSE_NO_DPD = get_NMSE(y0,x);
NMSE_ref_DPD = get_NMSE(y0_DPD,x);
%%
% plot_PM(x,y0,y_OPLS_DPD);
% plot_AM(x,y0,y_OPLS_DPD);
figure(2);
set(gcf, 'Color', 'w');
plot_AM_PM(x,y0,y_OPLS_DPD);
set(gca, 'Color', 'w');
