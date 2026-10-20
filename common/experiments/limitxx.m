function [PAPR,RMSin,b]=limitxx(x,a)

PAPR = 20*log10(max(abs(x))*sqrt(length(x))/norm(x)); 
RMSin = 10*log10( norm(x/a)^2/50/length(x/a)) + 30; %convert V-rms to dBm
b = PAPR+RMSin;
end