function verify_blind_theta_B()
% 验证：带宽 B（可用信息多少）对盲辨识 theta 精度的影响。
% 固定 K、SNR，扫 B，看 theta 误差如何变化。

clear; close all;

Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;

r_true=0.8; theta_true=0.6*pi;
b=[r_true^2,-2*r_true*cos(theta_true),1];
a=[1,-2*r_true*cos(theta_true),r_true^2];

f=(-N/2:N/2-1).'/N; w=2*pi*f;

% 循环谱先验（全频带算，mask 在拟合时再套）
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

th_grid=linspace(0,pi,3000);
K=64; snr=30;

% 在线估计 S_y^a（一次性算好，不随 B 变）
Sy=zeros(N,1);
for k=1:K
    s=gen(Nfft,Ncp,Nsym,M,1000+k);
    y=filter(b,a,s);
    n=(randn(N,1)+1j*randn(N,1))/sqrt(2);
    y=y+n*sqrt(mean(abs(y).^2))/sqrt(10^(snr/10));
    Y=fftshift(fft(y));
    Sy=Sy+Y.*conj(circshift(Y,dk));
end
Sy=Sy/K;
R_est=Sy./Ss;

B_list=[0.05,0.1,0.15,0.2,0.25,0.35,0.45,0.5];
err_B=zeros(size(B_list));

for bi=1:numel(B_list)
    B=B_list(bi);
    mask=abs(f)<=B/2;
    wgt=abs(Ss).*mask;
    score=zeros(size(th_grid));
    for ti=1:numel(th_grid)
        th=th_grid(ti);
        phi=allpass_phase(w,r_true,th);
        dphi=phi-circshift(phi,dk);
        R_model=exp(1j*dphi);
        score(ti)=real(wgt.' * (R_est.*conj(R_model)));
    end
    [~,imax]=max(score);
    theta_hat=th_grid(imax);
    err_B(bi)=abs(wrapToPi(theta_hat-theta_true));
    fprintf('B=%.2f (带内 %4d 点): theta_err=%.4f rad\n', B, nnz(mask), err_B(bi));
end

figure('Name','theta 误差 vs 带宽'); clf;
semilogy(B_list,err_B,'o-','LineWidth',1.5,'MarkerSize',8); grid on;
xlabel('带宽 B（归一化）'); ylabel('theta 辨识误差 (rad)');
title(sprintf('带宽对盲辨识 theta 精度的影响 (K=%d, SNR=%d dB)',K,snr));
end

function phi=allpass_phase(w,r,th)
    D=1-2*r*cos(th)*exp(-1j*w)+r^2*exp(-1j*2*w);
    phi=-2*w-2*angle(D);
end

function s=gen(Nfft,Ncp,Nsym,M,seed)
    rng(seed);
    data=randi([0 M-1],Nfft,Nsym);
    sym=qammod(data,M,'UnitAveragePower',true);
    td=ifft(sym,Nfft);
    s=zeros((Nfft+Ncp)*Nsym,1);
    for n=1:Nsym
        seg=[td(end-Ncp+1:end,n); td(:,n)];
        s((n-1)*(Nfft+Ncp)+1:n*(Nfft+Ncp))=seg;
    end
    s=s/sqrt(mean(abs(s).^2));
end
