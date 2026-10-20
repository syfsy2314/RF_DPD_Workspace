function D = geth_even_D(x,K,Q,h2,s)
%K = K+1;
% l = (Q-s+1)/s;
% d = zeros(length(x),s*(l+1));
% D = zeros(length(x),s);
% m = 1;
% for ll = 0:l
%     for n = (0+s*ll):(s-1+s*ll)
%       d(:,m) = x(:,1+K*n:K+K*n)* h2(1+K*ll:K+K*ll);
%       m = m+1;
%     end
% end
% 
% for a = 1:s
%     for ll = 0:l
%          D(:,a) = D(:,a) + d(:,ll*s+a);
%     end
% end
%%
% s为个数 l为组数
L = (Q-s+1)/s;
D = zeros(length(x),s);
d = zeros(length(x),(L+1)*s);
m = 1;
for l = 0:L
    for q = l*s:l*s+s-1
        for k =0:K-1
            d(:,m) = d(:,m) + h2((K)*l+k+1)*x(:,K*q+k+1);
        end
        m=m+1;
    end
end
%D有s列
for a = 1:s
    for l = 0 : L
        D(:,a) = D(:,a)+d(:,l*s+a);
    end
end
end