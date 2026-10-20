function [y_SBL,gain] = get_SBL(K,M,X,H_x,x,y_bf,N)
    n = K*(M+1);%系数个数
    % 初始化
    alpha = ones(n,1);
    sigma = var(y_bf)./10;
    % 迭代
    for i=1:N
        A = inv(diag(alpha));
        SIGMA =  inv(H_x' * H_x./sigma + inv(A)) ; 
        u = SIGMA * H_x' * y_bf./sigma; %#ok<MINV>
        alpha = 1./(abs(u).^2+diag(SIGMA));
        sigma = (norm(y_bf - H_x*u).^2 + trace(H_x' * H_x * SIGMA))/N;%#ok<MINV>
    end
    y_SBL = X * u;
    gain = rms(y_SBL)/rms(x);
end