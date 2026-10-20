function x = get_x
    GainBB=1; %%% GainBB between 0 and 1.
    [gBB,xsupply,PARAM] = IMS2017_generate_signal(GainBB);
    Fs=PARAM.CLK.CLKBB;
    fs = Fs /1e6;
    f_ofdm1 = PARAM.MOD.f_ofdm1;
    f_ofdm2 = PARAM.MOD.f_ofdm2;
    % nfft = 1024;
    N_samples=PARAM.L.LBB;
    t=(0:1:N_samples-1)';
    % 归一化
    PeakReal = max(abs(real(gBB)));  
    PeakImag = max(abs(imag(gBB)));
    Peak = max([PeakReal PeakImag]);
    x = double(gBB./Peak);
    [x1,x2]=demod_dual_BB(x,PARAM);
    P1 = abs(x1).^2;
    P2 = abs(x2).^2;
    Ps = P1+P2;
    %----CFR------------------------
    P_th = 0.3;
    [x_cfr, x_clip] = dual_band_clipper(x1,x2,P_th,PARAM);
    [x_cfr1,x_cfr2]=demod_dual_BB(x_cfr,PARAM);
    P1_cfr = abs(x_cfr1).^2;
    P2_cfr = abs(x_cfr2).^2;
    Ps_cfr = P1_cfr+P2_cfr;
    % 归一化
    PeakReal = max(abs(real(x_cfr2)));  
    PeakImag = max(abs(imag(x_cfr2)));
    Peak = max([PeakReal PeakImag]);
    % x_cfr2 = double(x_cfr2./Peak);
    x_ofdm =x_cfr2;
    x = x_ofdm + 0.0001*randn(20000,1);
    %x = gBB;
end