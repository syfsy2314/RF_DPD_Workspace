%% 测试自适应停止PLS
%%
clear;
close all;
fs = 1000;
nfft = 1024;
window=hann(nfft);
%% 第一个FIR滤波器
BW = 42;
%UNTITLED 返回离散时间滤波器对象。
Fs = 1000;  % Sampling Frequency
Fpass = BW+5;              % Passband Frequency
Fstop = BW+10;              % Stopband Frequency
% Dpass = 0.57501127785;  % Passband Ripple
Rp_dB = 0.5;        % 通带波纹，单位 3dB 大概是0.17
Dpass = (10^(Rp_dB/20) - 1) / (10^(Rp_dB/20) + 1);
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H1_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h1 = H1_FIR.Numerator;
% [EOF]
%% 第二个FIR滤波器
%UNTITLED 返回离散时间滤波器对象。
Fs = 1000;
Fpass = BW;
Fstop = BW+3;
Dpass = 0.0057501127785;
Dstop = 0.0001;
dens  = 20;
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
b  = firpm(N, Fo, Ao, W, {dens});
H2_FIR = dsp.FIRFilter( ...
    'Numerator', b);
h2 =  H2_FIR.Numerator;
%%
load C:\Users\admin\Desktop\测试\第四次\x.mat;
load C:\Users\admin\Desktop\测试\2026.5.28\test1\y0.mat;
gain_rmspa = rms(y0)/rms(x);
%%
% OMUX模拟滤波
y_f = H1_FIR(y0);
[y_f,x,~] = align_circular_signal(y_f,x);
% FIR滤波
y_bf = H2_FIR(y_f);
h = h2;
%% 初始化
%===========================================
% MP
K = 7;
M = 5;
L = K*(M+1);
H_x = geth_even_BLMP(x,K,M,h);%中频0 带限基函数
%===========================================
%% 提取线性分量并作投影
Len = 1;
x_lin = H_x(:,Len);
X_nl_raw = H_x(:,(Len+1):end);
% 2. 移走“大山”：对输出 y_bf 进行线性剥离
% 计算最佳线性系数 (最小二乘投影) 这一步算出的 w_lin 是物理上真实的线性增益
w_lin = (x_lin' * x_lin) \ (x_lin' * y_bf);
% 计算纯线性输出
y1 = x_lin * w_lin;
y_lin_fit = y1;
gain_rmsl = rms(y_lin_fit)/rms(x);
% 计算非线性残差 (这就是 PLS 要拟合的目标)
y_resid = y_bf - y_lin_fit;
gain_rms_res = rms(y_resid)/rms(x);
% 3. 清洗特征：对 X_nl_raw 进行去线性化
% 把每一列非线性特征里包含的线性成分减掉。
% 计算 X_nl_raw 在 x_lin 上的投影系数 (结果是 1x19 的行向量)
proj_coeffs = (x_lin' * x_lin) \ (x_lin' * X_nl_raw);
% 剥离线性成分
X_nl_clean = X_nl_raw - x_lin * proj_coeffs;
%% 进入PLS流程
%% ===================== 基本参数 =====================
T0 = size(H_x,2);              %% 【修改】避免后面 T 被潜变量矩阵覆盖
S = 0.5;
lamda = 1;
O = round(S * T0);             %% 原来固定 K=O，现在作为最大潜变量数
Kmax = O;                      %% 【修改】最大允许潜变量个数，而不是最终固定个数
X0 = X_nl_clean;
Y0 = y_resid;
%% ===================== 自适应停止参数 =====================
T_delta = 9e-4;                %% 【新增】边际残差下降阈值，可调：1e-4 ~ 1e-2
kappa_G_max = 1e3;               %% 【新增】潜变量矩阵条件数阈值，可调：1e3 ~ 1e6
A_min = 1;                     %% 【新增】至少保留的潜变量个数，防止过早停止
%% ===================== 数据标准化 =====================
N = size(X0,1);
x_mu = mean(X0,1);
x_sigma = std(X0,0,1);
x_sigma(x_sigma == 0) = eps;   %% 【新增】防止除零

y_mu = mean(Y0,1);
y_sigma = std(Y0,0,1);
if y_sigma == 0
    y_sigma = eps;             %% 【新增】防止除零
end

X = (X0 - x_mu) ./ x_sigma;
Y = (Y0 - y_mu) ./ y_sigma;

%% ===================== 初始化 =====================
W = zeros(size(X,2), Kmax);
P = zeros(size(X,2), Kmax);
Q = zeros(1, Kmax);

Z_score = zeros(N, Kmax);      %% 【新增】保存每个潜变量得分向量 t
Delta_list = nan(1, Kmax);     %% 【新增】保存每个潜变量的边际收益
Kappa_list = nan(1, Kmax);     %% 【新增】保存潜变量矩阵条件数

X_res = X;
Y_res = Y;

E_prev = real(Y_res' * Y_res); %% 【新增】当前残差能量
A = 0;                         %% 【新增】最终自适应得到的潜变量个数
stop_reason = 'Reach Kmax';    %% 【新增】记录停止原因

%% ===================== 自适应 PLS 潜变量提取 =====================
for k = 1:Kmax

    %% ---------- 1. 提取当前PLS方向 ----------
    w = X_res' * Y_res;
    w = w / norm(w);
    %% ---------- 2. 计算得分向量 ----------
    t = X_res * w;
    %% ---------- 3. 计算载荷 ----------
    p = X_res' * t / (t' * t);
    q = Y_res' * t / (t' * t);

    %% ---------- 4. 候选残差更新 ----------
    X_res_new = X_res - t * p';
    Y_res_new = Y_res - t * q;

    %% ---------- 5. 计算新增潜变量的边际收益 ----------
    E_new = real(Y_res_new' * Y_res_new);

    Delta_k = (E_prev - E_new) / (E_prev + eps);  %% 【新增】残差下降比例
    Delta_list(k) = Delta_k;

    %% ---------- 6. 计算加入当前潜变量后的病态程度 ----------
    % 临时加入当前候选 t；若条件数超限，当前 t 将在下一步被拒绝。
    Z_try = [Z_score(:,1:k-1), t];
    G_try = Z_try' * Z_try;

    kappa_G_k = cond(G_try);                %% 【修改】约束 H' * H 的条件数，更贴合后续求逆
    Kappa_list(k) = kappa_G_k;
    %% ---------- 7. 自适应停止准则 ----------
    % 【新增】如果新增潜变量收益太小，或者潜变量矩阵过于病态，则停止
    % 注意：前 A_min 个潜变量强制保留，避免一开始就停止
    if k > A_min
        if kappa_G_k >= kappa_G_max %% 留的少
            stop_reason = 'Normal-equation matrix is ill-conditioned';
            break;
        end
        if Delta_k <= T_delta %% 留的多
            stop_reason = 'Marginal residual gain is too small';
            break;
        end
    end

    %% ---------- 8. 保留当前潜变量 ----------
    W(:,k) = w;
    P(:,k) = p;
    Q(k) = q;
    Z_score(:,k) = t;

    X_res = X_res_new;
    Y_res = Y_res_new;
    E_prev = E_new;

    A = k;                     %% 【新增】当前保留的潜变量个数
end

%% ===================== 截断到最终自适应潜变量个数 =====================
if A == 0
    error('No latent variable is selected. Please reduce T_delta or check input data.');
end

W = W(:,1:A);
P = P(:,1:A);
Q = Q(:,1:A);
Z_score = Z_score(:,1:A);

Delta_list = Delta_list(1:A);
Kappa_list = Kappa_list(1:A);
fprintf('当前带宽为 = %dMHz\n', BW);
fprintf('自适应潜变量数 A = %d\n', A);
fprintf('Stopping reason: %s\n', stop_reason);

%% ===================== 计算PLS转换矩阵 =====================
% 原代码：R = W * inv(P' * W);
% 【修改】避免显式 inv，提高数值稳定性
R = W / (P' * W);
%% ===================== 降维 =====================
% 原代码：T = X * R;
% 【修改】用 Z 表示潜变量矩阵，避免和前面的 T0 混淆
Z = X * R;
%% ===================== 低维岭回归更新权重 =====================
I = eye(A);                    %% 【修改】维度由自适应 A 决定，不再是固定 K
% 原代码：C = inv(T' * T + lamda * I) * T' * Y;
% 【修改】避免 inv
C = (Z' * Z + lamda * I) \ (Z' * Y);
%% ===================== 还原到原始非线性BLMP系数 =====================
B_std = R * C;
% 还原尺度
B = B_std ./ (x_sigma.') * y_sigma;
%% 恢复 能不能剥多一点
% 计算 投影系数
X = geth_even_1d(x,K,M);
proj_coeff = (X(:,Len)' * X(:,Len)) \ (X(:,Len)' * X(:,(Len+1):end));
% 剥离线性成分
X_res = X(:,(Len+1):end) - X(:,Len) * proj_coeff;
% 恢复
y_OPLS = X_res*B +  X(:,Len)*w_lin;
gain_rmsOPLS = rms(y_OPLS)/rms(x);
%%
S = 0.25;
lamda = 1;
[y_LD_PLS_ref,gain_rms_LD_PLS_ref] = get_OPLS(K,M,X,x,h,y_bf,S,lamda);
NESE_LD_PLS_ref = get_NMSE(y_LD_PLS_ref,y0);
fprintf('最优潜变量数 A_ref = %d\n', round(S * 41));
%%
NESE_LD_PLS = get_NMSE(y_OPLS,y0);
