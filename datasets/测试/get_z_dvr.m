function [z] = get_z_dvr(K,M,x,pa_out)
    B = [0.1 0.3 0.5 0.7 0.9 1];
    H_y0 = geth_even_1d_DVR(pa_out/gain_peak,K,M,B);%中频0
    H_x = geth_even_1d_DVR(x,K,M,B);%中频0
    a0 = H_y0 \ x;
    z = H_x * a0;
end