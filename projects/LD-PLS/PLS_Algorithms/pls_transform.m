function B = pls_transform(X0,Y0,K) %返回系数 
    %标准化 新数据=（原数据- 均值）/ 标准差
    N = size(X0,1);
    X = (X0 - mean(X0))./std(X0);
    Y = (Y0 - mean(Y0))./std(Y0);
    W = zeros(size(X,2), K);
    P = zeros(size(X,2), K);
    Q = zeros(1, K);         % Y载荷向量（单输出）
    X_res = X;
    Y_res = Y;      % 初始化Y残差
    for k = 1:K
        w = X_res' * Y_res;%提取方向向量
        w = w / norm(w);%标准化
        t = X_res * w;%得分向量
        p = X_res' * t / (t'*t);%载荷矩阵 最小二乘
        q = Y_res' * t / (t' * t);  % Y载荷，单输出为标量
        X_res = X_res - t * p';%更新残差
        Y_res = Y_res - t * q;       % 更新Y残差
        W(:,k) = w;
        P(:,k) = p;
        Q(k) = q;
    end
    R = W * inv(P' * W);%转换矩阵R
    %降维
    T = X * R;
    % 降维系数 显然加上岭回归是有用的
    C0 = inv(T' * T ) * T' * Y;
    Yc = T * C0;
    % delta = (Yc - Y)'*(Yc - Y)/N;
    % lamda = K*delta;
    % lamda = 5;
    I = eye(K);
    C = inv(T' * T + lamda * I) * T' * Y;
    %C = inv(T' * T) * T' * Y;
    % 还原系数 
    B = R * C;
    % 还原尺度
    B = B./(std(X0)')*std(Y0);
end