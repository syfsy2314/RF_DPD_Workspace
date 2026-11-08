function [F,S]=F_matrix(x_ofdm,h_number)
N=length(x_ofdm);
F=zeros(length(x_ofdm),h_number);
S=zeros(h_number-1,h_number);
% if N > h_number
    for n=0:h_number-1 
        x=[zeros(n,1);x_ofdm(1:end-n)]; 
        F(:,n+1)=x;
    end
    for i=1:h_number-1
        y=[zeros(1,i),x_ofdm(end:-1:end-(h_number-i)+1).'];
        S(i,:)=y;
    end
% end