function y = run_mg(x,K,Q,h1,h2,s)
    % H=geth_even_1d(x,K,Q);
    % T = geth_even_T(H,K,Q,h1,s);
    % y = T*h2;
    H = geth_even_1d(x,K,Q);
    D = geth_even_D(H,K,Q,h2,s);
    y = D * h1;
end