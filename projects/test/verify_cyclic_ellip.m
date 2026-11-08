function verify_cyclic_ellip()
% 扩展验证：循环谱相位差方法对 IIR 椭圆滤波器（幅度+相位）是否依然成立。
% 核心：R = S_y^a / S_s^a = H(f)H*(f-a)，比较循环谱估计的 R_est 与真实 R_true。

clear;

Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;          % dk = 循环频率 a 的 bin 偏移

% IIR 椭圆滤波器（低通，通带内有幅度滚降 + 相位变化）
Wp=0.30; Ws=0.45; Rp=0.5; Rs=40;
[n,Wn]=ellipord(Wp,Ws,Rp,Rs);
[b,a]=ellip(n,Rp,Rs,Wn);
fprintf('椭圆阶数 n=%d\n', n);

f=(-N/2:N/2-1).'/N;

% 输入循环谱 S_s^a（已知 x，K0 次平均）
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

% 输出循环谱 S_y^a（K 次平均）
K=200; Sy=zeros(N,1);
for k=1:K
    s=gen(Nfft,Ncp,Nsym,M,1000+k);
    y=filter(b,a,s);
    Y=fftshift(fft(y));
    Sy=Sy+Y.*conj(circshift(Y,dk));
end
Sy=Sy/K;

R_est=Sy./Ss;                              % 循环谱比值（盲估计）

% 真实 R_true = H(f)H*(f-a)
H=fftshift(freqz(b,a,N,'whole'));
R_true=H.*conj(circshift(H,dk));

% 带内比较，且只在循环谱能量足够大的频点（权重 |Ss|）
mask = abs(f)<=0.30 & abs(f-dk/N)<=0.30;
wgt = abs(Ss) .* mask;
wgt = wgt / max(wgt);
good = wgt > 0.05;                          % 好频点（循环谱能量 > 5% 峰值）
in = find(good);

% 相位差诊断：真实相位差范围（是否超过 pi 造成缠绕）
dphi_true = angle(R_true(in));
fprintf('带内可用频点数: %d, 好频点数: %d\n', nnz(mask), numel(in));
fprintf('R_true 相位差范围: [%.0f, %.0f] deg\n', ...
    min(dphi_true)*180/pi, max(dphi_true)*180/pi);

e_phi_w = sum(wgt(in).*abs(wrapToPi(angle(R_est(in))-dphi_true))) / sum(wgt(in));
fprintf('R 相位差误差(加权): %.3f deg\n', e_phi_w*180/pi);

% 幅度只比较相对形状（各自归一化，消除整体增益）
a_est = abs(R_est(in)); a_true = abs(R_true(in));
a_est = a_est/mean(a_est); a_true = a_true/mean(a_true);
e_mag = mean(abs(a_est-a_true));
fprintf('R 幅度形状误差(归一化): %.4f\n', e_mag);
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
