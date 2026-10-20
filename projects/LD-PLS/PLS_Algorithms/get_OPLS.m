function [y_OPLS,gain_rms_all] = get_OPLS(K,M,X,x,h,y_bf,S,lamda)
H_x = geth_even_BLMP(x,K,M,h);
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
%%
T = size(H_x,2);
O = round (S * T);
B = pls_B(X_nl_clean,y_resid,O,lamda);
%% 恢复 能不能剥多一点
% 计算 投影系数 
proj_coeff = (X(:,Len)' * X(:,Len)) \ (X(:,Len)' * X(:,(Len+1):end));
% 剥离线性成分
X_res = X(:,(Len+1):end) - X(:,Len) * proj_coeff;
% 恢复
y_OPLS = X_res*B +  X(:,Len)*w_lin;
gain_rms_all = rms(y_OPLS)/rms(x);
end