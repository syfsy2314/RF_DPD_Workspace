function diag_ellip_phase()
% 相位准则辨识 Wn + 补偿验证
% (n,Rp,Rs) 已知，fminsearch 优化 Wn 使补偿后与 x 逐样本同相

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d); y0 = y0(:);

Wp = 0.10; Ws = 0.13; Rp = 0.5; Rs = 80;
[n, Wn_true] = ellipord(Wp, Ws, Rp, Rs);
[b_ell, a_ell] = ellip(n, Rp, Rs, Wn_true);

y0a = align_to(x, y0);
y_f = filter(b_ell, a_ell, y0a);

N = numel(x);
Yf = fft(y_f);
orth = @(z) rms(imag(z .* conj(x))) / rms(z);

fprintf('真值 Wn=%.4f, n=%d\n', Wn_true, n);

% fminsearch 优化 Wn（相位准则：正交分量最小）
cost = @(wn) orth(compensate(Yf, wn, n, Rp, Rs, N));
opts = optimset('Display','off','TolX',1e-7,'MaxIter',500);
wn0 = 0.08;
wn_opt = fminsearch(cost, wn0, opts);
fprintf('fminsearch 优化 Wn = %.4f (真值 %.4f)\n', wn_opt, Wn_true);

% 补偿（优化 Wn / 真值 Wn）
y_opt  = compensate(Yf, wn_opt,   n, Rp, Rs, N);
y_true = compensate(Yf, Wn_true,  n, Rp, Rs, N);
fprintf('优化Wn补偿: 正交分量 %.2f dB, 相对y0 NMSE %.2f dB\n', ...
    20*log10(orth(y_opt)), get_NMSE(y_opt, y0a));
fprintf('真值Wn补偿: 正交分量 %.2f dB, 相对y0 NMSE %.2f dB\n', ...
    20*log10(orth(y_true)), get_NMSE(y_true, y0a));
end

function y = compensate(Yf, Wn, n, Rp, Rs, N)
    [bt, at] = ellip(n, Rp, Rs, Wn);
    Ht = freqz(bt, at, N, 'whole');
    y = ifft(Yf .* exp(-1j*angle(Ht)));
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
