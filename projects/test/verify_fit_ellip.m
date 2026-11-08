function verify_fit_ellip()
% 验证：用非线性优化，从循环谱比值 R_est 参数化拟合椭圆滤波器的 b,a。
% 目标 max score = sum w(f) * real(R_est(f) * conj(H(f)H*(f-a)))
% 用真值附近初始化，先验证优化能否收敛回真值。

clear;
Nfft=32; Ncp=8; Nsym=50; M=4;
Ts=Nfft+Ncp; N=Nsym*Ts; dk=Nsym;
f=(-N/2:N/2-1).'/N;

% 椭圆真值
[b0,a0]=ellip(4,0.5,40,0.30);

% 输入循环谱 S_s^a
K0=400; Ss=zeros(N,1);
for k=1:K0
    s=gen(Nfft,Ncp,Nsym,M,k);
    S=fftshift(fft(s));
    Ss=Ss+S.*conj(circshift(S,dk));
end
Ss=Ss/K0;

% 输出循环谱 S_y^a
K=200; Sy=zeros(N,1);
for k=1:K
    s=gen(Nfft,Ncp,Nsym,M,1000+k);
    y=filter(b0,a0,s);
    Y=fftshift(fft(y));
    Sy=Sy+Y.*conj(circshift(Y,dk));
end
Sy=Sy/K;
R_est=Sy./Ss;

mask = abs(f)<=0.30 & abs(f-dk/N)<=0.30;
wgt = abs(Ss).*mask;

% 参数：x = [b(1:5); a(2:5)]，a(1)=1 固定
x_true = [b0(:); a0(2:end).'];

% 目标函数（嵌套，捕获 R_est,wgt,f,N,dk）
function val = negscore(x)
    b = x(1:numel(b0));
    a = [1; x(numel(b0)+1:end)];
    H = fftshift(freqz(b,a,N,'whole'));
    R_model = H.*conj(circshift(H,dk));
    val = -sum(wgt.*real(R_est.*conj(R_model)));
end

% 真值附近初始化 + 扰动
rng(7);
x_init = x_true + 0.05*randn(size(x_true));

opts = optimset('Display','off','MaxIter',2000,'MaxFunEvals',20000);
x_hat = fminsearch(@negscore, x_init, opts);

% 评估
b_hat = x_hat(1:numel(b0));
a_hat = [1; x_hat(numel(b0)+1:end)];
H_true = fftshift(freqz(b0,a0,N,'whole'));
H_hat  = fftshift(freqz(b_hat,a_hat,N,'whole'));
mask2 = abs(f)<=0.30;
% 归一化 DC 增益比较带内相位/幅度
c_true = H_true(find(mask2,1)); c_hat = H_hat(find(mask2,1));
H_true = H_true/c_true; H_hat = H_hat/c_hat;
e_phi = mean(abs(wrapToPi(angle(H_hat(mask2))-angle(H_true(mask2)))));
e_mag = mean(abs(abs(H_hat(mask2))-abs(H_true(mask2))));
fprintf('椭圆 b,a 拟合结果（带内）:\n');
fprintf('  相位误差: %.3f deg\n', e_phi*180/pi);
fprintf('  幅度误差: %.4f\n', e_mag);
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
