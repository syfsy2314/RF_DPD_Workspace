function nmse = get_NMSE(yk,yd)
    yd_nmse = yd/sqrt( mean( abs(yd).^2 ) );
    yk_nmse = yk/sqrt( mean( abs(yk).^2 ) );
    den = mean( abs(yd_nmse).^2 );
    num =  mean( abs(yk_nmse-yd_nmse).^2);
    nmse = 10*log10(num/den);
end