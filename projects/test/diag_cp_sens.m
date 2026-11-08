function diag_cp_sens()
% 诊断：CP 占比（循环平稳强度）对循环谱相位差估计精度的影响。
% 用全通滤波器（无幅度、小群时延），隔离 CP 占比这一个变量。

clear;
Nfft=256; Nsym=20; M=4;
r=0.8; th=0.6*pi;
b=[r^2,-2*r*cos(th),1]; a=[1,-2*r*cos(th),r^2];

cp_list=[8,16,32,64,128];
fprintf('CP  Ncp占比  相位差误差(deg)\n');
for ci=1:numel(cp_list)
    Ncp=cp_list(ci);
    Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;
    f=(-N/2:N/2-1).'/N;

    K0=300; Ss=zeros(N,1);
    for k=1:K0
        s=gen(Nfft,Ncp,Nsym,M,k);
        S=fftshift(fft(s));
        Ss=Ss+S.*conj(circshift(S,dk));
    end
    Ss=Ss/K0;

    K=300; Sy=zeros(N,1);
    for k=1:K
        s=gen(Nfft,Ncp,Nsym,M,1000+k);
        y=filter(b,a,s);
        Y=fftshift(fft(y));
        Sy=Sy+Y.*conj(circshift(Y,dk));
    end
    Sy=Sy/K;
    R_est=Sy./Ss;

    H=fftshift(freqz(b,a,N,'whole'));
    R_true=H.*conj(circshift(H,dk));

    mask=abs(f)<=0.4 & abs(f-dk/N)<=0.4;
    wgt=abs(Ss).*mask; wgt=wgt/max(wgt);
    in=find(wgt>0.05);
    e_phi=sum(wgt(in).*abs(wrapToPi(angle(R_est(in))-angle(R_true(in)))))/sum(wgt(in));
    fprintf('%2d   %.2f      %.3f\n', Ncp, Ncp/Ts, e_phi*180/pi);
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
