workspaceRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
datasetDir = fullfile(workspaceRoot, 'datasets');
testDir = fullfile(datasetDir, '测试', '2026.11.24', 'test2');
captureDir = fullfile(testDir, 'MP');
[~,~,x] = import_iq_binary(fullfile(datasetDir, 'data_38400_v0.6_test1.dat'), 38400);
gain = max(abs(x));
x = x./gain;
%% 导入采集的输出信号
[~,~,y0] = import_iq_binary(fullfile(captureDir, 'y0.bin'), 38400);
[y0,x,~] = align_circular_signal(y0,x);
y0 = y0./max(abs(y0));
% gain_rmspa = rms(y0)/rms(x);
% y0 = y0./gain_rmspa;
% gain = max(y0)/max(x);
%% 查看数据
nfft=1024;
fs = 1000;
FS = 1000e6;
N = 38400;
window=hann(nfft);

figure(1);
[p1, f1] = pwelch(x, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')
hold on;
[p5, f5] = pwelch(y0, window, nfft/2, nfft ,fs);
plot(f5 - 1/2 * fs, 10 * log10(fftshift(p5*100)), 'b-')
hold on;
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
grid on;
evm = get_EVM(y0,x);
%
figure(2);
t=1:38400;
plot(t,x);
hold on
plot(t,y0);
%% 导入预失真信号 转为可发送信号
data = load(fullfile(testDir, 'y0.mat'), 'y0');
y0 = data.y0;
z = get_z_mp(x,y0,7,5);
z = z * gain;
figure(1);
[p1, f1] = pwelch(z, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')
hold on;
%发送
SAVE_NAME = fullfile(captureDir, 'z0.dat');
fid = fopen(SAVE_NAME, 'wb');

for i = 1:length(z)
    fwrite(fid, real(z(i)), 'int16');
    fwrite(fid, imag(z(i)), 'int16');
end

fclose(fid);

%% bin转dat文件
folderPath = captureDir;
fileName = "y0_dpd";
fmat = fullfile(folderPath, fileName + ".mat");
fbin = fullfile(folderPath, fileName + ".bin");
[~,~,y0_dpd] = import_iq_binary(fbin,38400);
save(fmat, 'y0_dpd'); % 仅保存 myData 这一个变量
copyfile(fbin, fullfile(folderPath, fileName + ".dat"));
y = y0_dpd;
[y,x,~] = align_circular_signal(y,x);
get_EVM(y,x)
figure(1);
[p1, f1] = pwelch(x, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')
hold on;
[p5, f5] = pwelch(y0, window, nfft/2, nfft ,fs);
plot(f5 - 1/2 * fs, 10 * log10(fftshift(p5*100)), 'b-')
hold on;
gain1 = max(abs(y))/max(abs(x));
[p5, f5] = pwelch(y/gain1, window, nfft/2, nfft ,fs);
plot(f5 - 1/2 * fs, 10 * log10(fftshift(p5*100)), 'g-')
hold on;
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
grid on;
ACPR1=comm.ACPR('SampleRate',fs*1e6,...
    'MainChannelFrequency',0,...
    'MainMeasurementBandwidth',100e6,...
    'AdjacentMeasurementBandwidth',100e6,...
    'AdjacentChannelOffset',[-100e6,100e6],...
    'MainChannelPowerOutputPort',true,...
    'AdjacentChannelPowerOutputPort',true);
ACPR_DPD = step(ACPR1,y);