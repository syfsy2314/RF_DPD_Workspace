fs_dac=1000e6;
[x,Fs]=NR_2024_Signal_Generator(50e6,4096,15e3);
%star_symbol=NR_2024_Signal_Demodulator(x,60e6,4096,15e3);
calculate_PAPR(x);
% nfft=4096;
% window=hann(4096);
%[x,Fs]=NR_2024_Signal_Generator(100e6,4096,30e3);
%compare
figure;

I = real(x);
Q = imag(x);
IQ_1 = transpose(I + 1i * Q);
pwelch(IQ_1, [], [], [], 'centered',Fs);
figure;
%
change2dac=round(Fs/fs_dac);
dac_point=38400;


xx=transpose(x);
x=zeros(1,dac_point);

for i=1:dac_point
    x(i)=xx(i*change2dac);
end

x = x/max(abs(x));
x = 0.6*x;
%0.322 0.287 0.256 0.228   0.213    0.19  0.169  0.151 0.12
%                         -21.8dbm  -22.3 -22.8  -23.3 -24.3
x = round(x*(2^(16-1)-10));%放大到16位
calculate_PAPR(x);
SAVE_NAME = "D:\研一科研\NR信号发生器及毫米波PA模型\data_38400_v0.6_50MHz_test1.dat";
fid = fopen(SAVE_NAME, 'wb');

for i = 1:length(x)
    fwrite(fid, real(x(i)), 'int16');
    fwrite(fid, imag(x(i)), 'int16');
end

fclose(fid);
I = real(x);
Q = imag(x);

IQ = I + 1i * Q;

pwelch(IQ, [], [], [], 'centered',fs_dac);
% 
% [p1, f1] = pwelch(IQ, window, nfft/2, nfft ,fs_dac);
% plot(f1 - 1/2 * fs_dac, 10 * log10(fftshift(p1)), 'r-')%红色
% hold on;


%%
% x2 =x;
% x2=x2*exp(1i*pi/3);
% %x2 = round(x2*(2^(16-1)-10));%放大到16位
% 
% SAVE_NAME = "E:\zynq_test_data\rfdc_dac\test_1t1r_1GHz\MIMO\data_phase_20000_v0.1_test1.dat";
% fid = fopen(SAVE_NAME, 'wb');
% 
% for i = 1:length(x2)
%     fwrite(fid, real(x2(i)), 'int16');
%     fwrite(fid, imag(x2(i)), 'int16');
% end
% 
% fclose(fid);
% figure;
% I = real(x);
% Q = imag(x);
% 
% IQ = I + 1i * Q;
% 
% pwelch(IQ, [], [], [], 'centered',fs_dac);
function papr_value = calculate_PAPR(x)
% 计算信号的峰均功率比 (PAPR)
% 输入参数：
%   x: 输入信号（可以是复数）
% 输出：
%   papr_value: PAPR值（单位：dB）

% 计算瞬时功率
instant_power = abs(x).^2;

% 计算峰值功率
peak_power = max(instant_power);

% 计算平均功率
avg_power = mean(instant_power);

% 计算PAPR（线性值）
papr_linear = peak_power / avg_power;

% 转换为dB
papr_value = 10 * log10(papr_linear);

% 显示结果
fprintf('峰值功率: %.6f\n', peak_power);
fprintf('平均功率: %.6f\n', avg_power);
fprintf('PAPR: %.2f dB\n', papr_value);
end
%%
function recovered_constellation = NR_2024_Signal_Demodulator(signal, bandwidth, FFT, ScSpacing)
%% 参数设置（与生成器一致）
OSF = 16;                                            % 过采样因子
FFT_size = FFT * OSF;                               % IFFT/FFT大小
BW_utilization = 0.95;                              % 占用的信道带宽比例
NumSc = ceil(bandwidth * BW_utilization / ScSpacing) + 1; % 活跃子载波数
if mod(NumSc, 2)
    NumSc = NumSc - 1;
end
Fs = FFT_size * ScSpacing;                          % 采样率
nOFDMSymbols = 7 * 2;                                % OFDM符号数
M_QAM = 64;                                         % 子载波调制阶数
CP_duration = 144 * OSF;                            % CP长度
Nt = 16 * OSF;                                      % 窗过渡长度

%% 信号分割与去窗处理
symbol_length = CP_duration + FFT_size;             % 每个符号总长度
recovered_symbols = zeros(FFT_size, nOFDMSymbols);

% 计算符号起始位置
symbol_starts = zeros(1, nOFDMSymbols);
symbol_starts(1) = 1;
for i = 2:nOFDMSymbols
    % 每个符号的起始位置（考虑重叠部分）
    symbol_starts(i) = symbol_starts(i-1) + symbol_length - Nt;
end

% 提取并处理每个符号
for i = 1:nOFDMSymbols
    % 提取当前符号（包括CP和窗过渡）
    if i < nOFDMSymbols
        symbol = signal(symbol_starts(i):symbol_starts(i) + symbol_length - 1);
    else
        % 最后一个符号包含额外的后缀
        symbol = signal(symbol_starts(i):symbol_starts(i) + symbol_length + Nt - 1);
    end
    
    % 去窗处理
    if i > 1
        % 去除前一个符号的后缀影响
        symbol(1:Nt) = symbol(1:Nt) - signal(symbol_starts(i-1) + symbol_length:...
                                             symbol_starts(i-1) + symbol_length + Nt - 1);
    end
    
    % 去除循环前缀
    symbol = symbol(CP_duration+1:end);
    
    % 如果是最后一个符号，去掉额外的后缀
    if i == nOFDMSymbols
        symbol = symbol(1:FFT_size);
    end
    
    recovered_symbols(:, i) = symbol;
end

%% FFT变换
fft_output = fft(recovered_symbols) / sqrt(FFT_size);

%% 提取活跃子载波
Sc_Index = [1:NumSc/2, FFT_size-NumSc/2+1:FFT_size];
recovered_qam_symbols = fft_output(Sc_Index, :);

%% 归一化处理
% 由于发送端使用了单位功率归一化，这里进行相应处理
recovered_qam_symbols = recovered_qam_symbols / sqrt(mean(abs(recovered_qam_symbols(:)).^2));

%% 绘制星座图
figure;
scatter(real(recovered_qam_symbols(:)), imag(recovered_qam_symbols(:)), 10, 'filled');
xlabel('实部 (I)');
ylabel('虚部 (Q)');
title('恢复的星座图 (64-QAM)');
grid on;
axis equal;
hold on;

% 绘制理想星座点
ideal_constellation = qammod(0:M_QAM-1, M_QAM, 'UnitAveragePower', true);
scatter(real(ideal_constellation), imag(ideal_constellation), 100, 'r', 'x', 'LineWidth', 1.5);
legend('恢复的符号', '理想星座点', 'Location', 'best');

% 计算并显示EVM
evm = sqrt(mean(abs(recovered_qam_symbols(:) - ideal_constellation(randi([1, M_QAM], numel(recovered_qam_symbols), 1))).^2));
fprintf('估计的EVM: %.2f%%\n', evm * 100);

%% 返回恢复的符号
recovered_constellation = recovered_qam_symbols(:);
end
function recovered_constellation = NR_2024_Signal_Demodulator_fixed(signal, bandwidth, FFT, ScSpacing)
%% 参数设置
OSF = 16;
FFT_size = FFT * OSF;
BW_utilization = 0.95;
NumSc = ceil(bandwidth * BW_utilization / ScSpacing) + 1;
if mod(NumSc, 2)
    NumSc = NumSc - 1;
end
nOFDMSymbols = 7 * 2;
M_QAM = 64;
CP_duration = 144 * OSF;
Nt = 16 * OSF;

%% 计算正确的符号起始位置
% 每个符号的基础长度（CP + FFT_size）
symbol_basic_length = CP_duration + FFT_size;  % 2304 + 16384 = 18688

% 第一个符号的起始和结束位置
symbol_starts = zeros(1, nOFDMSymbols);
symbol_ends = zeros(1, nOFDMSymbols);

% 第一个符号：有后缀
symbol_starts(1) = 1;
symbol_ends(1) = symbol_basic_length + Nt;  % 18688 + 256 = 18944

% 计算后续符号的起始位置
for i = 2:nOFDMSymbols
    % 后续符号与前一个符号有Nt个样点的重叠
    symbol_starts(i) = symbol_starts(i-1) + symbol_basic_length - Nt;
    if i < nOFDMSymbols
        symbol_ends(i) = symbol_starts(i) + symbol_basic_length - 1;  % 普通符号
    else
        symbol_ends(i) = symbol_starts(i) + symbol_basic_length + Nt - 1;  % 最后一个符号
    end
end

% 验证最后一个索引是否超出信号长度
if symbol_ends(end) > length(signal)
    error('计算错误：符号索引超出信号长度。信号长度：%d，需要：%d', ...
        length(signal), symbol_ends(end));
end

%% 提取和解调每个符号
recovered_symbols = zeros(FFT_size, nOFDMSymbols);

for i = 1:nOFDMSymbols
    % 提取符号
    if i == 1
        % 第一个符号：有后缀
        symbol = signal(symbol_starts(i):symbol_ends(i));
        symbol = symbol(1:symbol_basic_length);  % 去掉后缀
    elseif i < nOFDMSymbols
        % 中间符号：有前后重叠
        symbol = signal(symbol_starts(i):symbol_ends(i));
    else
        % 最后一个符号：有前缀
        symbol = signal(symbol_starts(i):symbol_ends(i));
    end
    
    % 去窗处理
    if i > 1
        % 前一个符号的后缀
        prev_symbol_end = symbol_starts(i) - 1;
        prev_postfix = signal(prev_symbol_end - Nt + 1:prev_symbol_end);
        
        % 从当前符号的前Nt个点中减去前一个符号的后缀
        symbol(1:Nt) = symbol(1:Nt) - prev_postfix;
    end
    
    % 移除CP
    symbol = symbol(CP_duration + 1:end);
    
    % 如果是最后一个符号，去掉额外的后缀
    if i == nOFDMSymbols
        symbol = symbol(1:FFT_size);
    end
    
    % 存储
    recovered_symbols(:, i) = symbol(1:FFT_size);
end

%% FFT和解调
Sc_Index = [1:NumSc/2, FFT_size - NumSc/2 + 1:FFT_size];
recovered_qam = zeros(NumSc, nOFDMSymbols);

for i = 1:nOFDMSymbols
    fft_output = fft(recovered_symbols(:, i)) / sqrt(FFT_size);
    recovered_qam(:, i) = fft_output(Sc_Index);
end

%% 绘制星座图
recovered_qam_flat = recovered_qam(:) / sqrt(mean(abs(recovered_qam(:)).^2));

figure;
scatter(real(recovered_qam_flat), imag(recovered_qam_flat), 20, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
xlabel('实部 (I)');
ylabel('虚部 (Q)');
title(sprintf('恢复的星座图 (共%d个点)', length(recovered_qam_flat)));
grid on;
axis equal;

% 添加理想星座点
ideal_constellation = qammod(0:M_QAM-1, M_QAM, 'UnitAveragePower', true);
hold on;
scatter(real(ideal_constellation), imag(ideal_constellation), 100, 'r', 'x', 'LineWidth', 1.5);
legend('恢复的符号', '理想星座点', 'Location', 'best');

% 计算EVM
evm = calculate_evm(recovered_qam_flat, M_QAM);
fprintf('EVM: %.2f%%\n', evm * 100);
fprintf('平均功率: %.4f\n', mean(abs(recovered_qam_flat).^2));

recovered_constellation = recovered_qam_flat;
end

function evm = calculate_evm(received_symbols, M)
% 计算EVM
ideal_symbols = qammod(0:M-1, M, 'UnitAveragePower', true);
% 找到每个接收符号对应的最近理想星座点
[~, indices] = min(abs(received_symbols - ideal_symbols.'), [], 2);
decoded = ideal_symbols(indices).';
evm = sqrt(mean(abs(received_symbols - decoded).^2)) / sqrt(mean(abs(decoded).^2));
end