function verify_blind_rtheta()
% 完全盲辨识：同时辨识全通滤波器参数 (r, theta)，二维网格搜索。
%
% 全通滤波器相位 phi(w;r,th) 依赖两个参数，从循环谱相位差
%   dphi = phi(f;r,th) - phi(f-a;r,th)
% 用加权相关在 (r,th) 二维网格上最大化。

clear; close all;

Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;

r_true=0.8; theta_true=0.6*pi;
b=[r_true^2,-2*r_true*cos(theta_true),1];
a=[1,-2*r_true*cos(theta_true),r_true^2];

f=(-N/2:N/2-1).'/N; w=2*pi*f;

% 循环谱先验 S_s^a
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

mask=abs(f)<=0.25/2;
wgt=abs(Ss).*mask;

% 二维参数网格
r_grid=linspace(0.4,0.98,60);
th_grid=linspace(0,pi,400);

SNR_list=[40,30,25,20,15,10];
K_list=[64,256];
err=zeros(numel(SNR_list),numel(K_list));

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
        R_est=Sy./Ss;

        % 二维网格搜索
        score=zeros(numel(r_grid),numel(th_grid));
        for ri=1:numel(r_grid)
            r=r_grid(ri);
            D=1-2*r*cos(th_grid).*exp(-1j*w)+r^2*exp(-1j*2*w);  % N x T
            phi=-2*w-2*angle(D);                                % N x T
            dphi=phi-circshift(phi,dk);
            R_model=exp(1j*dphi);
            score(ri,:)=real(wgt.' * (R_est.*conj(R_model)));
        end
        [~,idx]=max(score(:));
        [ri,ti]=ind2sub(size(score),idx);
        r_hat=r_grid(ri); theta_hat=th_grid(ti);
        err(si,ki)=abs(wrapToPi(theta_hat-theta_true));
        fprintf('SNR=%2d dB, K=%3d : r=%.3f(%d), theta=%.4f(%d), dtheta=%.4f rad\n', ...
            snr, K, r_hat, r_hat==r_true, theta_hat, theta_hat==theta_true, err(si,ki));
    end
end
fprintf('r_true=%.3f, theta_true=%.4f rad\n', r_true, theta_true);
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
