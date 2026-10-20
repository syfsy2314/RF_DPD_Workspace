function Lb = get_Lb(H,Kb,nb)
K=size(H,2);
Index=2:(Kb+1);
for i = 1:K
    for j = 1:K
        d(i,j) = dist(H(:,i).' , H(:,j));  % 计算欧几里得距离
    end
end
[B,I]=sort(d,2);
B = B(:,Index);%保留2：kb列
I = I(:,Index);%保留2：kb列
B = B.';
I = I.';
R = repmat((1 : K),[Kb, 1]);%它的作用是将行向量 1:K 垂直复制 Kb 次，生成一个 Kb × K 的矩阵。
% R(:) 是一种常用的语法，用于将矩阵 R 转换为列向量。它的作用是将矩阵的所有元素按列优先顺序重新排列成一个列向量。
A = sparse(R(:), I(:), B(:), K, K);
A = A + A.';
% A = sparse(R(:), I(:), B(:), K, K) 是 MATLAB 中用于 构造稀疏矩阵 的一种高效方式。它的作用是基于给定的 行索引 R、列索引 I 和值 B，生成一个大小为 K×K 的稀疏矩阵 A。
%上述步骤找出所有欧式距离中最小的两个值，以及索引，把其他元素清零
tb = nb.*(sum(A,"all")./nnz(A));
Wb = exp(-A./tb);
D = sum(Wb,2);
Db = diag(D);
Lb = Db - Wb;
end