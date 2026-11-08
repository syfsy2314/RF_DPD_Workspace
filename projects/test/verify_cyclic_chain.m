function verify_cyclic_chain()
% 正确设定：窄带 OFDM + 带限观测，循环谱只恢复带内相位（相位差链），不拟合系数。
% 核心：R = S_y^a/S_s^a = H(f)H*(f-a)，相位差 D(f)=angle(H(f))-angle(H(f-a))
%       相位差链：angle(H(f)) = angle(H(f0)) + 累加 D，只在带内累加。

clear;
Nfft=512; Ncp=16; Nsym=8; M=128;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;          % dk = (1/Ts)*N = Nsym
f=(-N/2:N/2-1).'/N;

% 信号带宽（窄带，只占 ~0.1 归一化）
B_s = 0.05;

% 椭圆滤波器：通带覆盖信号，阻带在带外
[b0,a0]=ellip(6,0.5,50,0.15);

% 输入循环谱 S_s^a（alpha=符号率）
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,B_s,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

% 输出循环谱 S_y^a（带限观测 y_f = filter(H, s)）
K=400; Sy=zeros(N,1);
for k=1:K
    s=gen(Nfft,Ncp,Nsym,M,B_s,1000+k);
    y=filter(b0,a0,s);
    Y=fftshift(fft(y));
    Sy=Sy+Y.*conj(circshift(Y,dk));
end
Sy=Sy/K;

R_est=Sy./Ss;

% 真实滤波器相位
H_true=fftshift(freqz(b0,a0,N,'whole'));
phi_true=unwrap(angle(H_true));

% 带内可辨识区间：f 与 f-a 都落在信号带内
a_freq = dk/N;
band = B_s + a_freq;                        % 信号带 + 循环频率裕量
mask = abs(f)<=B_s & abs(f-a_freq)<=B_s;
wgt = abs(Ss).*mask; wgt=wgt/max(wgt);
good = wgt>0.05; in=find(good);

% 相位差 D(f) = angle(R_est)
D = angle(R_est);

% 相位差链：从带内中心 f0=0 向两侧累加（步长 a_freq）
% 直接对 D 在 f 网格上做累加重建（近似）
[~,i0] = min(abs(f));
% 向右（正频率）累加
phi_rec = zeros(N,1);
idx = i0:dk:N/2;                             % 沿正频率，步长 dk
for m=1:numel(idx)-1
    phi_rec(idx(m+1)) = phi_rec(idx(m)) + D(idx(m));
end
% 向左（负频率）累加
idx = i0:-dk:1;
for m=1:numel(idx)-1
    phi_rec(idx(m+1)) = phi_rec(idx(m)) - D(idx(m+1));
end

% 比较带内重建相位 vs 真实相位（去常数）
in_band = abs(f)<=B_s;
phi_rec_ub = phi_rec(in_band);
phi_tru_ub = phi_true(in_band);
phi_rec_ub = phi_rec_ub - mean(phi_rec_ub);
phi_tru_ub = phi_tru_ub - mean(phi_tru_ub);
e_phi = mean(abs(phi_rec_ub - phi_tru_ub));

fprintf('窄带 B_s=%.3f, 符号率 a=%.5f, 带内频点数=%d\n', B_s, a_freq, nnz(in_band));
fprintf('带内相位重建误差(去常数): %.3f deg\n', e_phi*180/pi);
end

function s=gen(Nfft,Ncp,Nsym,M,B_s,seed)
    rng(seed);
    Nused = max(2, round(Nfft*B_s*2));       % 占用子载波数（双边带宽）
    Nused = 2*floor(Nused/2);
    data = randi([0 M-1], Nused, Nsym);
    sym = qammod(data, M, 'UnitAveragePower', true);
    td = zeros(Nfft, Nsym);
    idx = [1:Nused/2, Nfft-Nused/2+1:Nfft];   % 中间子载波
    td(idx, :) = sym;
    td = ifft(td, Nfft);
    s = zeros((Nfft+Ncp)*Nsym, 1);
    for n=1:Nsym
        seg=[td(end-Ncp+1:end,n); td(:,n)];
        s((n-1)*(Nfft+Ncp)+1:n*(Nfft+Ncp))=seg;
    end
    s=s/sqrt(mean(abs(s).^2));
end
