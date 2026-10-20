function T = geth_even_T(x,K,Q,h1,s)
%K = K+1;
% l = (Q-s+1)/s;
% T = zeros(length(x),K*(l+1));
% m = 1;
% for ll = 0:l
%     for q = ll*s:ll*s+s-1
%         for k = 1:K
%            tk(:,m) = x(:,k+K*q) * h1(q-ll*s+1);
%            m = m+1;
%         end
%         T(:,1+ll*K:K+ll*K) =  T(:,1+ll*K:K+ll*K) + tk; 
%         m = 1;
%     end
% end
%%
L =  (Q+1-s)/s;
T = zeros(length(x),K*(L+1));
m=1;
for l = 0:L
    for k=0:K-1
        for q = l*s:l*s+s-1
            T(:,m)=T(:,m)+h1(q-l*s+1)*x(:,q*K+k+1);
        end
        m=m+1;
    end
end

end