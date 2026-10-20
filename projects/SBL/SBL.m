clear;clc;close all;
nfft=1024;
fs = 200;
FS = 200e6;
N = 20000;
window=hanning(nfft); 
%% 生成OFDM
% GainBB=1; %%% GainBB between 0 and 1.
% [gBB,xsupply,PARAM] = IMS2017_generate_signal(GainBB);
% Fs=PARAM.CLK.CLKBB;
% fs = Fs /1e6;
% f_ofdm1 = PARAM.MOD.f_ofdm1;
% f_ofdm2 = PARAM.MOD.f_ofdm2;
% nfft = 1024;
% N_samples=PARAM.L.LBB;
% t=[0:1:N_samples-1]';
% % 归一化
% PeakReal = max(abs(real(gBB)));  
% PeakImag = max(abs(imag(gBB)));
% Peak = max([PeakReal PeakImag]);
% x = double(gBB./Peak);
% [x1,x2]=demod_dual_BB(x,PARAM);
% P1 = abs(x1).^2;
% P2 = abs(x2).^2;
% Ps = P1+P2;
% %----CFR------------------------
% P_th = 0.3;
% [x_cfr, x_clip] = dual_band_clipper(x1,x2,P_th,PARAM);
% [x_cfr1,x_cfr2]=demod_dual_BB(x_cfr,PARAM);
% P1_cfr = abs(x_cfr1).^2;
% P2_cfr = abs(x_cfr2).^2;
% Ps_cfr = P1_cfr+P2_cfr;
% % 归一化
% PeakReal = max(abs(real(x_cfr2)));  
% PeakImag = max(abs(imag(x_cfr2)));
% Peak = max([PeakReal PeakImag]);
% % x_cfr2 = double(x_cfr2./Peak);
% x_ofdm =x_cfr2;
% % x = x_ofdm;
% x = x_ofdm + 0.0001*randn(20000,1);
load x.mat
%% 过功放
% PeakReal = max(abs(real(x)));  
% PeakImag = max(abs(imag(x)));
% Peak = max([PeakReal PeakImag]);
% PAPR = 20*log10(max(abs(x))*sqrt(length(x))/norm(x));
% % 过网上功放
% [PAPR_o,RMSin_o,b]=limitxx(x,15);
% [pa_out, RMSout] = RFWebLab_PA_meas_v1_2(x, RMSin_o);%无DPD输出
% [pa_out,x,~] = align_circular_signal(pa_out,x);
% M = 3;%记忆深度
% K  = 5;%奇数阶的个数 ，代表的是 1, 3, 5阶
% load pa_coeff
% pacoe=20;pdcoe=1;    %功放放大倍数pacoe,预失真器放大倍数pdcoe;
% pa_out= pacoe * volterra(x, pa_coeff, 3, 2); %功放的输出 %the desired PA output
load pa_out.mat
y0 = pa_out;
gain_rmspa = rms(y0)/rms(x);%单PA增益
%% 滤波器 
Hd = LPF_8;%低通滤波 带宽8MHz
h = Hd.numerator;
h = h';
R = H_matrix(h,N);
y_bf = R * y0;
%初始化
K = 5;
M = 3;
H_x = geth_even_BLMP(x,K,M,h);%中频0 带限基函数
X = geth_even_1d(x,K,M);
a = inv(ctranspose(H_x) * H_x) * ctranspose(H_x) * y_bf;
y_BLMP = H_x * a;
gain_rmsBLMP = rms(y_BLMP)/rms(x);
%% SBL
n = K*(M+1);%系数个数
% 初始化
alpha = ones(n,1);
sigma = var(y_bf)./10;
% 迭代
for i=1:5
    A = inv(diag(alpha));
    SIGMA =  inv(H_x' * H_x./sigma + inv(A)) ; 
    u = SIGMA * H_x' * y_bf./sigma; %#ok<MINV>
    for j=1:n
        alpha = 1./(abs(u).^2+SIGMA(j,j));
    end
    sigma = (norm(y_bf - H_x*u).^2 + trace(H_x' * H_x * SIGMA))/N;%#ok<MINV>
end
y_SBL = X * u;
% gain_rmsSBL = rms(y_SBL)/rms(x);
gain_peakSBL = abs(max(y_SBL))/abs(max(x));
H_SBL = geth_even_1d(y_SBL/gain_peakSBL,K,M);
A1 = pinv(H_SBL) * x;
Z1 = X * A1;
% 过网上功放
[PAPR_o,RMSin_o,b]=limitxx(Z1,15);
[y1, RMSout] = RFWebLab_PA_meas_v1_2(Z1, RMSin_o);%无DPD输出
[y1,x,~] = align_circular_signal(y1,x);
gain_rmsSBL = rms(y1)/rms(x);
%% MG-BLMP
K = 5; Q = 3; s = 2;
% 初始化
h2 = zeros(K*((Q-s+1)/s+1),1);
h2(1)=1;
D = geth_even_D(H_x,K,Q,h2,s);
h1 = pinv(D) * y_bf;
T = geth_even_T(H_x,K,Q,h1,s);
h2 = pinv(T) * y_bf;
%第二遍
D = geth_even_D(H_x,K,Q,h2,s);
h1 = pinv(D) * y_bf;
T = geth_even_T(H_x,K,Q,h1,s);
h2 = pinv(T) * y_bf;
y_MG = run_mg(x,K,Q,h1,h2,s);
% gain_rmsMG = rms(y_MG)/rms(x);
gain_peakMG = abs(max(y_MG))/abs(max(x));
H_MG = geth_even_1d(y_MG/gain_peakMG,K,M);
A2 = pinv(H_MG) * x;
Z2 = X * A2;
% 过网上功放
[PAPR_o,RMSin_o,b]=limitxx(Z2,15);
[y2, RMSout] = RFWebLab_PA_meas_v1_2(Z2, RMSin_o);%无DPD输出
[y2,x,~] = align_circular_signal(y2,x);
gain_rmsMG = rms(y2)/rms(x);
%% RR
delta = (R*y_BLMP - y_bf)'*(R*y_BLMP - y_bf)/N;
L = K*(M+1);
u = L*delta/N;
u = 0.00015;
a3 = inv(ctranspose(H_x) * H_x + u * eye(L)) * ctranspose(H_x) * y_bf;
y_RR = X * a3;
% gain_rmsRR = rms(y_RR)/rms(x);
gain_peakRR = abs(max(y_RR))/abs(max(x));
H_RR = geth_even_1d(y_RR/gain_peakRR,K,M);
A5 = pinv(H_RR) * x;
Z5 = X * A5;
% 过网上功放
[PAPR_o,RMSin_o,b]=limitxx(Z5,15);
[y5, RMSout] = RFWebLab_PA_meas_v1_2(Z5, RMSin_o);%无DPD输出
[y5,x,~] = align_circular_signal(y5,x);
gain_rmsRR = rms(y5)/rms(x);
%% MR
Kb = 2;
nb = 0.2;
a1 = 0.00001;
Lb = get_Lb(H_x,Kb,nb);
%
Kg = 20;
ng = 0.2;
a2 = 0.00001;
g = 1000;
G = H_x(1:g,:);
Lg = get_Lg(G,Kg,ng);
%
pa4= inv(ctranspose(H_x) * H_x + a1 * Lb + a2 * G' * Lg * G) * ctranspose(H_x) * y_bf;
y_MR = X * pa4;
% gain_rmsMR = rms(y_MR)/rms(x);
gain_peakMR = abs(max(y_MR))/abs(max(x));
H_MR = geth_even_1d(y_MR/gain_peakMR,K,M);
A4 = pinv(H_MR) * x;
Z4 = X * A4;
% 过网上功放
[PAPR_o,RMSin_o,b]=limitxx(Z4,15);
[y4, RMSout] = RFWebLab_PA_meas_v1_2(Z4, RMSin_o);%无DPD输出
[y4,x,~] = align_circular_signal(y4,x);
gain_rmsMR = rms(y4)/rms(x);
%% AMR
a = 5;
b = 0.0001;
B = [0.1 0.3 0.5 0.7 0.9 1];
Xa = geth_even_1d_DVR(x,5,3,B);
H_Xa = R * Xa;
I2 = eye(size(Xa,2));
Xb = geth_even_1d_DDR(x,7,2);
Tb = ctranspose(Xa)*Xa-ctranspose(Xa)*Xb*inv(ctranspose(Xb)*Xb)*ctranspose(Xb)*Xa;
pa5 = inv(ctranspose(H_Xa) * H_Xa + a * Tb + b * I2) * ctranspose(H_Xa) * y_bf;
y_AMR = Xa * pa5;
% gain_rmsAMR = rms(y_AMR)/rms(x);
gain_peakAMR = abs(max(y_AMR))/abs(max(x));
H_AMR = geth_even_1d(y_AMR/gain_peakAMR,K,M);
A3 = pinv(H_AMR) * x;
Z3 = X * A3;
% 过网上功放
[PAPR_o,RMSin_o,b]=limitxx(Z3,15);
[y3, RMSout] = RFWebLab_PA_meas_v1_2(Z3, RMSin_o);%无DPD输出
[y3,x,~] = align_circular_signal(y3,x);
gain_rmsAMR = rms(y3)/rms(x);
%% PCA
n = 10;
% 带限输入信号降维
data = zscore(H_x'*H_x);%数据标准化
% data = U'*U;
r = corrcoef(data);%相关系数矩阵，L*L维度 L里面找N个
% 利用相关系数矩阵进行主成分分析，vec1的第一列为r的第一特征向量
[vec1,lamda,b] = pcacov(r);%lamda为r的特征值
P = vec1(:,1:n);% n 为降维后的系数
U4 = H_x * P;
a4 = inv(ctranspose(U4) * U4 ) * ctranspose(U4) * y_bf;
% 原始输入信号降维
% data = zscore(H'*H);%数据标准化
y_PCA = X * P * a4;
% gain_rmsPCA = rms(y_PCA)/rms(x);
gain_peakPCA = abs(max(y_PCA))/abs(max(x));
H_PCA = geth_even_1d(y_PCA/gain_peakPCA,K,M);
A6 = pinv(H_PCA) * x;
Z6 = X * A6;
% 过网上功放
[PAPR_o,RMSin_o,b]=limitxx(Z6,15);
[y6, RMSout] = RFWebLab_PA_meas_v1_2(Z6, RMSin_o);%无DPD输出
[y6,x,~] = align_circular_signal(y6,x);
gain_rmsPCA = rms(y6)/rms(x);
%%
figure(1);
[p1, f1] = pwelch(y0/gain_rmspa, window, nfft/2, nfft ,fs);
plot(f1 - 1/2 * fs, 10 * log10(fftshift(p1*100)), 'r-')%红色
hold on;
[p4, f4] = pwelch(y1/gain_rmsSBL, window, nfft/2, nfft ,fs);
plot(f4 - 1/2 * fs, 10 * log10(fftshift(p4*100)), 'c-')%绿色
hold on;
[p6, f6] = pwelch(y2/gain_rmsMG, window, nfft/2, nfft ,fs);
plot(f6 - 1/2 * fs, 10 * log10(fftshift(p6*100)), 'b-')%绿色
hold on;
[p3, f3] = pwelch(y3/gain_rmsAMR, window, nfft/2, nfft ,fs);
plot(f3 - 1/2 * fs, 10 * log10(fftshift(p3*100)), 'g-')%绿色
hold on;
[p2, f2] = pwelch(y4/gain_rmsMR, window, nfft/2, nfft ,fs);
plot(f2 - 1/2 * fs, 10 * log10(fftshift(p2*100)), 'k-')%绿色
hold on;
[p5, f5] = pwelch(y5/gain_rmsRR, window, nfft/2, nfft ,fs);
plot(f5 - 1/2 * fs, 10 * log10(fftshift(p5*100)), 'm-')%绿色
hold on;
[p7, f7] = pwelch(y6/gain_rmsPCA, window, nfft/2, nfft ,fs);
plot(f7 - 1/2 * fs, 10 * log10(fftshift(p7*100)),'Color',[0.8,0.6,0.1]);
hold on;
xlabel('频率偏移 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
legend('无DPD','SBL','MG-BLMP','ExMR','AMR','RR','PCA');
grid on;
axis([-35 35 -60 10]);
%%
y_SBL = y1;
y_MG = y2;
y_AMR = y3;
y_MR = y4;
y_RR = y5;
y_PCA = y6;
y0_nmse = y0/sqrt( mean( abs(y0).^2 ) );
y1_nmse = y_SBL/sqrt( mean( abs(y_SBL).^2 ) );
y2_nmse = y_MG/sqrt( mean( abs(y_MG).^2 ) );
y3_nmse = y_AMR/sqrt( mean( abs(y_AMR).^2 ) );
y4_nmse = y_MR/sqrt( mean( abs(y_MR).^2 ) );
y5_nmse = y_RR/sqrt( mean( abs(y_RR).^2 ) );
y6_nmse = y_PCA/sqrt( mean( abs(y_PCA).^2 ) );
num1=mean( abs(y0_nmse-y1_nmse).^2);
num2=mean( abs(y0_nmse-y2_nmse).^2);
num3=mean( abs(y0_nmse-y3_nmse).^2);
num4=mean( abs(y0_nmse-y4_nmse).^2);
num5=mean( abs(y0_nmse-y5_nmse).^2);
num6=mean( abs(y0_nmse-y6_nmse).^2);
den=mean( abs(y0_nmse).^2 );
NMSE_SBL = 10*log10(num1/den);
NMSE_MG_BLMP = 10*log10(num2/den);
NMSE_AMR = 10*log10(num3/den);
NMSE_ExMR = 10*log10(num4/den);
NMSE_RR = 10*log10(num5/den);
NMSE_PCA = 10*log10(num6/den);
%% ACPR
ACPR1=comm.ACPR('SampleRate',fs*1e6,...
    'MainChannelFrequency',0,...
    'MainMeasurementBandwidth',10e6,...
    'AdjacentMeasurementBandwidth',8e6,...
    'AdjacentChannelOffset',[-10e6,10e6],...
    'MainChannelPowerOutputPort',true,...
    'AdjacentChannelPowerOutputPort',true);
ACPR_No_DPD = step(ACPR1,y0);
ACPR_SBL_DPD = step(ACPR1,y1);
ACPR_MG_DPD = step(ACPR1,y2);
ACPR_AMR_DPD = step(ACPR1,y3);
ACPR_MR_DPD = step(ACPR1,y4);
ACPR_RR_DPD = step(ACPR1,y5);
ACPR_PCA_DPD = step(ACPR1,y6);
