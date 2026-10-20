function [y] = volterra(x, coeff, order, memory)

%-- coeff is the coefficient of volterra
%-- order is  the odd order of volterra 
%-- memory is the memory depth of volterra

len_x = length(x);
y = zeros(len_x, 1);%
beta = zeros(len_x, memory + 1);
% len_coeff = length(coeff);

for num = (memory + 1): 1 : len_x
    for num_coeff = 1 : 1 : order * (memory + 1)
        memory_temp = floor((num_coeff - 1) / (order));                %��9��ϵ��Ϊ����memory_temp = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        orders_temp = mod(num_coeff - 1, order);                               %��9��ϵ��Ϊ��: orders_temp = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        
       y_temp(num, num_coeff) = coeff(num_coeff) * x(num - memory_temp)... 
                                 * abs(x(num - memory_temp))^orders_temp^2;   % -----ϵ�� * x(n-q) * |x(n-q)|^����
%         beta(num, memory_temp + 1) = beta(num, memory_temp + 1) + y_temp(num, num_coeff);     %��ǰ��������µĽ������ۼӣ�
        y(num, 1) = y(num, 1) + y_temp(num, num_coeff);                     %-----����ۼӹ���
    end    
end


