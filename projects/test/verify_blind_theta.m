function verify_blind_theta()
% 盲辨识全通滤波器参数 theta（r 已知），用循环谱相位差做参数化拟合。
%
% 核心：R(ω)=S_y^a(ω)/S_s^a(ω)=H(ω)H*(ω-a)=exp(j[phi(ω;θ)-phi(ω-a;θ)])
%       盲辨识 θ = 网格搜索最大化加权相关  Σ w(ω)·R_est(ω)·conj(R_model(ω;θ))
% 相比逐点重建相位曲线，参数化拟合用所有频点投票估计单个参数，更稳健。

clear; close all;

Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;       % dk = (1/Ts)*N

% 真值
r_true=0.8; theta_true=0.6*pi;
b=[r_true^2,-2*r_true*cos(theta_true),1];
a=[1,-2*r_true*cos(theta_true),r_true^2];

f=(-N/2:N/2-1).'/N;                    % fftshift 频率轴
w=2*pi*f;

% 循环谱先验 S_s^a(f)（数值平均，代表已知 OFDM 格式）
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

% 带限 + 权重（循环谱幅度大的频点权重高，近零点自然压低）
mask=abs(f)<=0.25/2;
wgt=abs(Ss).*mask;

% 参数网格
th_grid=linspace(0,pi,3000);

SNR_list=[40,30,25,20,15,10];
K_list=[16,64,256];
theta_err=zeros(numel(SNR_list),numel(K_list));

for si=1:numel(SNR_list)
    snr=SNR_list(si);
    for ki=1:numel(K_list)
        K=K_list(ki);
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
        R_est=Sy./Ss;                    % 复循环谱比值

        % 网格搜索 theta：最大化加权相关
        score=zeros(size(th_grid));
        for ti=1:numel(th_grid)
            th=th_grid(ti);
            phi=allpass_phase(w,r_true,th);
            dphi=phi-circshift(phi,dk);  % phi(f)-phi(f-a)
            R_model=exp(1j*dphi);
            score(ti)=real(wgt.' * (R_est.*conj(R_model)));
        end
        [~,imax]=max(score);
        theta_hat=th_grid(imax);
        theta_err(si,ki)=abs(wrapToPi(theta_hat-theta_true));
        fprintf('SNR=%2d dB, K=%3d : theta_hat=%.4f rad, err=%.4f rad\n', ...
            snr, K, theta_hat, theta_err(si,ki));
    end
end

figure('Name','盲辨识 theta 误差'); clf;
leg=cell(numel(K_list),1);
for ki=1:numel(K_list)
    semilogy(SNR_list,theta_err(:,ki),'o-','LineWidth',1.5,'MarkerSize',8); hold on;
    leg{ki}=sprintf('K=%d',K_list(ki));
end
xlabel('SNR (dB)'); ylabel('theta 辨识误差 (rad)'); grid on;
legend(leg,'Location','best');
title('盲辨识全通滤波器参数 theta');
end

function phi=allpass_phase(w,r,th)
    % 二阶全通滤波器相位响应
    % H(e^{jw}) = e^{-j2w} * conj(D)/D, D = 1 - 2r cos(th) e^{-jw} + r^2 e^{-j2w}
    % phi = -2w - 2*angle(D)
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
