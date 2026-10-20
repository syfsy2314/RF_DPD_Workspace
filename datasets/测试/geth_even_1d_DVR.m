function H=geth_even_1d_DVR(x,K,M,B)
s = (M+1) + (M+1) * K * 3;
H=zeros(length(x),s);
n=1;
for i = 0:M   
    x1 = [zeros(i,1);x(1:end-i)];
    H(:,n) = x1;
    n = n+1;
end
for i = 0:M   
    x1 = [zeros(i,1);x(1:end-i)];
    for k = 1:K   
            H(:,n) = abs(abs(x1)-B(k)).* exp(j * angle(x1));
            n = n+1;
    end
end
for i = 0:M   
    x1 = [zeros(i,1);x(1:end-i)];
    for k = 1:K   
            H(:,n) = abs(abs(x1)-B(k)).* exp(j * angle(x1)).* abs(x);
            n = n+1;
    end
end
for i = 0:M   
    x1 = [zeros(i,1);x(1:end-i)];
    for k = 1:K   
            H(:,n) = abs(abs(x1)-B(k)).*x;
            n = n+1;
    end
end
for i = 0:M   
    x1 = [zeros(i,1);x(1:end-i)];
    for k = 1:K   
            H(:,n) = abs(abs(x1)-B(k)).*x1;
            n = n+1;
    end
end
for i = 0:M   
    x1 = [zeros(i,1);x(1:end-i)];
    for k = 1:K   
            H(:,n) = abs(abs(x)-B(k)).*x1;
            n = n+1;
    end
end
end

