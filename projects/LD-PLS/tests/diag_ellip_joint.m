function diag_ellip_joint()
% 椭圆结构先验 + PA 物理结构先验 联合辨识诊断
% 阶段1/2：真实数据 AM/PM 量化 + 椭圆带限滤波器参数确定 + 生成 y_f

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d); y0 = y0(:);
K = 7; M = 5;

% 对齐（整数延迟 + 常数相位，不含频率相关群时延）
y0a = align_to(x, y0);

% ---- 阶段1：量化 AM/PM 失配（同相基函数覆盖 y0 的残差）----
H = geth_even_1d(x, K, M);
a_pa = (H'*H) \ (H'*y0a);
nmse_ampm = get_NMSE(H*a_pa, y0a);
fprintf('--- 阶段1: AM/PM 失配量化 ---\n');
fprintf('对齐后同相拟合 NMSE(=AM/PM+阶数残差): %.2f dB\n', nmse_ampm);
fprintf('PA 增益 rms(y0)/rms(x): %.4f\n', rms(y0a)/rms(x));

% ---- 阶段2：椭圆带限滤波器（统一 OMUX+覆盖）----
% 信号带宽 ±0.048(归一化)，通带取 0.10 留余量，阻带 0.13
Wp = 0.10; Ws = 0.13; Rp = 0.5; Rs = 80;
[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b_ell, a_ell] = ellip(n, Rp, Rs, Wn);
fprintf('--- 阶段2: 椭圆带限滤波器 ---\n');
fprintf('椭圆阶数 n=%d, 归一化截止 Wn=%.4f\n', n, Wn);

y_f = filter(b_ell, a_ell, y0a);
y_f = align_to(x, y_f);
fprintf('y_f 相对 y0 的 NMSE(带限损失): %.2f dB\n', get_NMSE(y_f, y0a));

% ---- 阶段3：可辨识性分析（已知 H 时 PA 系数的可辨识性）----
fprintf('--- 阶段3: 可辨识性分析 ---\n');
X = geth_even_1d(x, K, M);          % N x 42, 全带基函数
Hx = filter(b_ell, a_ell, X);        % 每列过椭圆 = H_theta X
rho = sum(abs(Hx).^2, 1) ./ sum(abs(X).^2, 1);   % 各列带内能量占比
rho_k = reshape(rho, K, M+1);        % 行=k(非线性阶), 列=m+1(延迟)
fprintf('各非线性阶 k 的带内能量占比 rho_k（7个值=延迟m=0..5）:\n');
for k = 1:K
    fprintf('  k=%d: %s\n', k, sprintf('%.3f  ', rho_k(k,:)));
end
s_hx = svd(Hx); s_hx = s_hx / s_hx(1);
s_x  = svd(X);  s_x  = s_x  / s_x(1);
fprintf('Hx 奇异值(带限) 前8=[%s]\n', sprintf('%.2g ', s_hx(1:8)));
fprintf('X  奇异值(全带) 前8=[%s]\n', sprintf('%.2g ', s_x(1:8)));
fprintf('Hx 条件数=%.3g,  X 条件数=%.3g\n', s_hx(1)/s_hx(end), s_x(1)/s_x(end));
fprintf('Hx >1e-3: %d/%d,   X >1e-3: %d/%d\n', ...
    nnz(s_hx > 1e-3), numel(s_hx), nnz(s_x > 1e-3), numel(s_x));

% ---- 阶段4a：已知 H 的覆盖+外推 误差地板 ----
fprintf('--- 阶段4a: 已知 H 的覆盖+外推误差地板 ---\n');
[U, S, V] = svd(Hx, 'econ');
sv = diag(S);
for r = [3 5 8 10 13 20 30 42]
    Ur = U(:,1:r); Sr = S(1:r,1:r); Vr = V(:,1:r);
    a_hat = Vr * (Sr \ (Ur' * y_f));
    y_hat = X * a_hat;
    fprintf('  覆盖+外推 r=%2d: NMSE = %.2f dB\n', r, get_NMSE(y_hat, y0a));
end
a_dir = (X'*X) \ (X'*y_f);
fprintf('  直接拟合(全带X->y_f): NMSE = %.2f dB\n', get_NMSE(X*a_dir, y0a));

% ---- 对照：正交化基函数（消除共线性），看带限奇异值谱是否仍有悬崖 ----
fprintf('--- 对照: 正交化基函数后带限可辨识性 ---\n');
[Qmat, ~] = qr(X, 0);
Hq = filter(b_ell, a_ell, Qmat);
sq = svd(Hq); sq = sq / sq(1);
fprintf('正交基 Hq >1e-3: %d/%d, 条件数=%.3g, 前8=[%s]\n', ...
    nnz(sq > 1e-3), numel(sq), sq(1)/sq(end), sprintf('%.2g ', sq(1:8)));

% ---- 正交化覆盖+外推 误差地板 ----
fprintf('--- 正交化覆盖+外推 误差地板 ---\n');
[Uq, Sq, Vq] = svd(Hq, 'econ');
for r = [5 10 20 30 42]
    c = Vq(:,1:r) * (Sq(1:r,1:r) \ (Uq(:,1:r)' * y_f));
    fprintf('  正交化覆盖+外推 r=%2d: NMSE = %.2f dB\n', r, get_NMSE(Qmat*c, y0a));
end
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
