function [y_PLS,gain] = get_PLS(K,M,X,x,h,y_bf,S)
H_x = geth_even_BLMP(x,K,M,h);
T = K * (M+1);
O = round (S * T);
B = pls_transform(H_x,y_bf,O);
y_PLS = X * B;
gain = rms(y_PLS)/rms(x);
end