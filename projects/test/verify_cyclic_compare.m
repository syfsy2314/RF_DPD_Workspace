function verify_cyclic_compare()
% 对照实验：全通 vs 椭圆，直接比较循环谱比值 R_est 与真实 R_true。
% 目的：定位椭圆失败是"实现 bug"还是"滤波器特性（群时延/幅度）"。

clear;
Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;
f=(-N/2:N/2-1).'/N;

% 输入循环谱 S_s^a
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

fprintf('--- 滤波器类型对照 ---\n');
for typ = 1:2
    if typ==1
        r=0.8; th=0.6*pi;
        b=[r^2,-2*r*cos(th),1]; a=[1,-2*r*cos(th),r^2];
        name='全通 r=0.8';
        band=0.50;
    else
        [~,Wn]=ellipord(0.30,0.45,0.5,40);
        [b,a]=ellip(4,0.5,40,Wn); name='椭圆 n=4';
        band=0.30;
    end

    % 输出循环谱 S_y^a
    K=200; Sy=zeros(N,1);
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

    mask = abs(f)<=band & abs(f-dk/N)<=band;
    wgt = abs(Ss).*mask; wgt=wgt/max(wgt);
    good = wgt>0.05; in=find(good);

    dphi_true = angle(R_true(in));
    e_phi = sum(wgt(in).*abs(wrapToPi(angle(R_est(in))-dphi_true)))/sum(wgt(in));
    fprintf('%s: 好频点=%d, R_true相位差范围=[%.0f,%.0f]deg, 相位误差=%.3f deg\n', ...
        name, numel(in), min(dphi_true)*180/pi, max(dphi_true)*180/pi, e_phi*180/pi);
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
