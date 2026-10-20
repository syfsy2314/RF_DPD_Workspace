function [y,gain] = get_MG(K,M,s,x,h,y_bf)
    H_x = geth_even_BLMP(x,K,M,h);%中频0 带限基函数
    % MG
    Q = M;
    % 初始化
    h2 = zeros(K*((Q-s+1)/s+1),1);
    h2(1)=1;
    D = geth_even_D(H_x,K,Q,h2,s);
    h1 = pinv(D) * y_bf;
    T = geth_even_T(H_x,K,Q,h1,s);
    h2 = pinv(T) * y_bf;
    %第二遍
    D = geth_even_D(H_x,K,Q,h2,s);
    h1 = pinv(D) * y_bf;
    T = geth_even_T(H_x,K,Q,h1,s);
    h2 = pinv(T) * y_bf;
    y = run_mg(x,K,Q,h1,h2,s);
    gain = rms(y)/rms(x);
end