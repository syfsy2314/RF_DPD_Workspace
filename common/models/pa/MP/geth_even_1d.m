function H=geth_even_1d(x,K,M)
s = K * (M+1);
H=zeros(length(x),s);
n=1;
for m=0:M
    x1=[zeros(m,1);x(1:end-m)];
    for k=1:K
            H(:,n)=(abs(x1).^(k-1)).*x1;
            n = n+1;
    end
end

end