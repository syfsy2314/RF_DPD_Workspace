function verify_ellip_split()
% 分步辨识验证：幅度用功率谱(alpha=0)，相位用循环谱(alpha=符号率)。
% |H(f)|^2 = S_y^0/S_s^0（普通功率谱比，无需循环平稳，应很可靠）
% 相位差  = angle(S_y^a/S_s^a)（循环谱，重点）

clear;
Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;
f=(-N/2:N/2-1).'/N;

[b0,a0]=ellip(4,0.5,40,0.30);

% 功率谱 S^0（alpha=0，普通功率谱）
K0=200; Ps=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ps=Ps+abs(S).^2;
end
Ps=Ps/K0;

K=200; Py=zeros(N,1);
for k=1:K
    s=gen(Nfft,Ncp,Nsym,M,1000+k);
    y=filter(b0,a0,s);
    Y=fftshift(fft(y));
    Py=Py+abs(Y).^2;
end
Py=Py/K;

% 幅度估计 |H| = sqrt(Py/Ps)
mask = abs(f)<=0.30;
H_true = fftshift(freqz(b0,a0,N,'whole'));
mag_est = sqrt(Py./Ps);
mag_est = mag_est/mag_est(find(mask,1));
mag_true = abs(H_true)/abs(H_true(find(mask,1)));
e_mag = mean(abs(mag_est(mask)-mag_true(mask))./mag_true(mask));
fprintf('幅度估计相对误差(功率谱比): %.4f\n', e_mag);

% 循环谱 S^a（alpha=符号率）
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

K=400; Sy=zeros(N,1);
for k=1:K
    s=gen(Nfft,Ncp,Nsym,M,1000+k);
    y=filter(b0,a0,s);
    Y=fftshift(fft(y));
    Sy=Sy+Y.*conj(circshift(Y,dk));
end
Sy=Sy/K;

R_est=Sy./Ss;
R_true=H_true.*conj(circshift(H_true,dk));

mask2 = abs(f)<=0.30 & abs(f-dk/N)<=0.30;
wgt = abs(Ss).*mask2; wgt=wgt/max(wgt);
in = find(wgt>0.05);
e_phi = sum(wgt(in).*abs(wrapToPi(angle(R_est(in))-angle(R_true(in)))))/sum(wgt(in));
fprintf('相位差估计误差(循环谱, K=400): %.3f deg\n', e_phi*180/pi);
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
