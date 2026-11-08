function diag_joint_id_diag()
% 诊断过参数化 + 交叉关系分离失败原因。

d = load('projects/LD-PLS/data/x.mat');  x = struct2array(d); x = x(:);
d = load('projects/LD-PLS/data/y0.mat'); y0 = struct2array(d); y0 = y0(:);
N = numel(x);

fs = 1000;
Wp = 35/(fs/2); Ws = 55/(fs/2); Rp = 0.5; Rs = 50;
[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b0, a0] = ellip(n, Rp, Rs, Wn);

y_f = align_to(x, filter(b0, a0, y0));

K = 7; M = 5; P = 60;
L = P - 1 + M; Ma = M + 1; Ld = L + 1;

% --- 1) 过参数化回归拟合质量 ---
Psi = zeros(N, Ld*K);
col = 0;
for l = 0:L
    xl = [zeros(l,1); x(1:end-l)];
    axl = abs(xl);
    for k = 1:K
        col = col + 1;
        Psi(:, col) = xl .* axl.^(k-1);
    end
end
dvec = (Psi'*Psi) \ (Psi'*y_f);
fprintf('过参数化回归拟合 y_f: NMSE = %.2f dB\n', get_NMSE(Psi*dvec, y_f));

% --- 2) d_1 与 h_true 的相似度 ---
D = reshape(dvec, Ld, K);
h_true = impz(b0, a0, P);
d1 = D(1:P, 1);
corr = abs(d1'*h_true) / (norm(d1)*norm(h_true));
fprintf('d_1 前%d点 与 h_true 相关: %.4f\n', P, corr);

% --- 3) A_mat 奇异值谱（判断解唯一性）---
nb_rows = (Ma + Ld - 1) * (K-1);
A_mat = zeros(nb_rows, K*Ma);
for k = 2:K
    r0 = (k-2)*(Ma + Ld - 1);
    T1 = convmat(D(:,1), Ma);
    Tk = convmat(D(:,k), Ma);
    A_mat(r0+1:r0+Ma+Ld-1, 1:Ma)        = -Tk;
    A_mat(r0+1:r0+Ma+Ld-1, (k-1)*Ma+1:k*Ma) = T1;
end
s = svd(A_mat);
fprintf('A_mat 奇异值谱（前5 / 末8）:\n');
fprintf('  前5: '); fprintf('%.2e ', s(1:5)); fprintf('\n');
fprintf('  末8: '); fprintf('%.2e ', s(end-7:end)); fprintf('\n');
fprintf('  条件数: %.2e\n', s(1)/s(end));
end

function T = convmat(d, len)
Ld = length(d);
T = zeros(Ld + len - 1, len);
for i = 1:len
    T(i:i+Ld-1, i) = d;
end
end

function ya = align_to(x, y)
    Hf = fft(x) .* conj(fft(y));
    [~, ind] = max(abs(ifft(Hf)));
    ya = circshift(y, -(ind-1));
    c = (x' * ya) / (x' * x);
    ya = ya * conj(c) / abs(c);
end
