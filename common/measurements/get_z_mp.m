function [z] = get_z_mp(x,y,K,M)
    X = geth_even_1d(x,K,M);
    gain_peak = abs(max(y))/abs(max(x));
    lambda = 1e-4;
    H = geth_even_1d(y/gain_peak,K,M);
    % A = (H' * H + lambda * eye(size(H,2))) \ (H' * x);
    % A = pinv(H) * x;
    A = H \ x;
    z = X * A;
end