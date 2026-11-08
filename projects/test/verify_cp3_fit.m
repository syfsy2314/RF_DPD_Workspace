function verify_cp3_fit()
% 验证：真实 CP 占比(3%)下，参数化拟合（单参数 theta，r 已知）能否聚合噪声。
% 若 CP 3% 下参数化拟合误差仍小，说明"参数化拟合"能救回低循环平稳性。

clear;
Nfft=512; Ncp=16; Nsym=8; M=128;      % 真实信号参数，CP 占比 3%
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;
f=(-N/2:N/2-1).'/N;

r_true=0.8; theta_true=0.6*pi;
b=[r_true^2,-2*r_true*cos(theta_true),1];
a=[1,-2*r_true*cos(theta_true),r_true^2];

% 输入循环谱
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

mask=abs(f)<=0.4 & abs(f-dk/N)<=0.4;
wgt=abs(Ss).*mask;

th_grid=linspace(0,pi,3000);
w=2*pi*f;

for K=[100,400,1000]
    Sy=zeros(N,1);
    for k=1:K
        s=gen(Nfft,Ncp,Nsym,M,1000+k);
        y=filter(b,a,s);
        Y=fftshift(fft(y));
        Sy=Sy+Y.*conj(circshift(Y,dk));
    end
    Sy=Sy/K;
    R_est=Sy./Ss;

    score=zeros(size(th_grid));
    for ti=1:numel(th_grid)
        th=th_grid(ti);
        D=1-2*r_true*cos(th)*exp(-1j*w)+r_true^2*exp(-1j*2*w);
        phi=-2*w-2*angle(D);
        dphi=phi-circshift(phi,dk);
        score(ti)=real(wgt.' * (R_est.*conj(exp(1j*dphi))));
    end
    [~,imax]=max(score);
    theta_hat=th_grid(imax);
    fprintf('K=%4d: theta_hat=%.4f, 误差=%.4f rad (%.2f deg)\n', ...
        K, theta_hat, abs(wrapToPi(theta_hat-theta_true)), abs(wrapToPi(theta_hat-theta_true))*180/pi);
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
