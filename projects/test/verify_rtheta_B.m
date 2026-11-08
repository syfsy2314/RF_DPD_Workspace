function verify_rtheta_B()
% 两参数 (r,theta) 辨识：扫带宽 B，看是否带宽增大后能解耦 r 和 theta。
% 固定 SNR=30, K=256，二维网格搜索。

Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;

r_true=0.8; theta_true=0.6*pi;
b=[r_true^2,-2*r_true*cos(theta_true),1];
a=[1,-2*r_true*cos(theta_true),r_true^2];

f=(-N/2:N/2-1).'/N; w=2*pi*f;

K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

% 有噪多段平均 Sy（不随 B 变）
K=256; snr=30; Sy=zeros(N,1);
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

r_grid=linspace(0.4,0.98,80);
th_grid=linspace(0,pi,400);

B_list=[0.15,0.2,0.25,0.3,0.35,0.4,0.5];
for bi=1:numel(B_list)
    B=B_list(bi);
    mask=abs(f)<=B/2;
    wgt=abs(Ss).*mask;
    score=zeros(numel(r_grid),numel(th_grid));
    for ri=1:numel(r_grid)
        r=r_grid(ri);
        D=1-2*r*cos(th_grid).*exp(-1j*w)+r^2*exp(-1j*2*w);
        phi=-2*w-2*angle(D);
        Rm=exp(1j*(phi-circshift(phi,dk)));
        score(ri,:)=real(wgt.'*(R_est.*conj(Rm)));
    end
    [~,idx]=max(score(:));
    [ri,ti]=ind2sub(size(score),idx);
    r_hat=r_grid(ri); theta_hat=th_grid(ti);
    fprintf('B=%.2f : r_hat=%.3f (真0.80), theta_hat=%.4f (真1.885), 误差=%.3f/%.3f\n', ...
        B, r_hat, theta_hat, abs(r_hat-r_true), abs(wrapToPi(theta_hat-theta_true)));
end
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
