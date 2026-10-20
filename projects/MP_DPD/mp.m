clear;clc;close all;
nfft=1024;
fs = 200;
FS = 200e6;
N = 20000;
window=hanning(nfft);
x = get_x;
%% 过功放
PeakReal = max(abs(real(x)));
PeakImag = max(abs(imag(x)));
Peak = max([PeakReal PeakImag]);
PAPR = 20*log10(max(abs(x))*sqrt(length(x))/norm(x));
% pa_out = get_paout_online(x,x);
M = 1;%记忆深度
K  = 7;%奇数阶的个数 ，代表的是 1, 3, 5阶
load pa_coeff.mat;
pacoe=20;pdcoe=1;    %功放放大倍数pacoe,预失真器放大倍数pdcoe;
pa_out= pacoe * volterra(x, pa_coeff, 3, 2); %功放的输出 %the desired PA output
gain_rmspa = rms(pa_out)/rms(x);%单PA增益
gain_peak = abs(max(pa_out))/abs(max(x));
%% MP DPD 间接学习结构 需要在线更新输出
% I = 8;%迭代次数
% H_x = geth_even_1d(x,K,M);%中频0
% z0 = x;%初始化
% A0 = zeros(K*(M+1),1);
% A0(1,1)=1;%初始化DPD系数
% yy=zeros(length(x),I);
% gain_peaky0 = abs(max(y0))/abs(max(x));
% H_y0 = geth_even_1d(y0/gain_peaky0,K,M);
% gain_rmsy0=zeros(1,I);
% for i=1:I
%     A0 = inv(H_y0'*H_y0)*(H_y0')*z0;%DPD系数%最小化x=H_y0*A0(DPD输出)的误差，只迭代了一次，所以u=x
%     z0 = H_x * A0;
%     yy(:,i) = pacoe * volterra(z0, pa_coeff, 3, 2);
%     gain_peaky0 = abs(max(yy(:,i)))/abs(max(x));
%     gain_rmsy0(1,i) = rms(yy(:,i))/rms(x);
%     H_y0 = geth_even_1d(yy(:,i)/gain_peaky0,K,M);
% % yy = pacoe * volterra(z0, pa_coeff, 3, 2);
% % [PAPR_o,RMSin_o,b]=limitxx(z0,15);
% % [pd_pa_out, RMSout] = RFWebLab_PA_meas_v1_2(z0, RMSin_o);%无DPD输出
% % [pd_pa_out,x,~] = align_circular_signal(pd_pa_out,x);
% % yy = pd_pa_out;
% % gain_rmsy0 = rms(yy)/rms(x);
% end
%% MP DPD 直接学习结构 加权 LS 算法迭代估计系数
H_y0 = geth_even_1d(pa_out/gain_peak,K,M);%中频0
H_x = geth_even_1d(x,K,M);%中频0
a0 = pinv(H_y0)*x;
z0 = H_x * a0;
%y0 = get_paout_online(z0,x);
y0= pacoe * volterra(z0, pa_coeff, 3, 2); %功放的输出 %the desired PA output
gain_rms = rms(y0)/rms(x);%单PA增益
%%
figure(1);
[p1, f1] = pwelch(x, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1)), 'r-')
hold on;
[p2, f2] = pwelch(pa_out/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f2 - 1/2 * fs, 10 * log10(fftshift(p2)), 'g-')
hold on;
[p3, f3] = pwelch(y0/gain_rms, window, nfft/2, nfft ,fs);
plot(f3 - 1/2 * fs, 10 * log10(fftshift(p3)), 'b-')
hold on;
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
legend('输入信号','无DPD','MP');
grid on;
%%
ACPR1=comm.ACPR('SampleRate',fs*1e6,...
    'MainChannelFrequency',0,...
    'MainMeasurementBandwidth',30e6,...
    'AdjacentMeasurementBandwidth',25e6,...
    'AdjacentChannelOffset',[-30e6,30e6],...
    'MainChannelPowerOutputPort',true,...
    'AdjacentChannelPowerOutputPort',true);
ACPR_No_DPD = step(ACPR1,pa_out);
ACPR_MP_DPD = step(ACPR1,y0);
%% nmse
NMSE_MP = get_NMSE(y0,x);
NMSE_NODPD = get_NMSE(pa_out,x);
%% EVM
EVM_MP = get_EVM(y0,x);
EVM_NODPD = get_EVM(pa_out,x);
%% 星座图
