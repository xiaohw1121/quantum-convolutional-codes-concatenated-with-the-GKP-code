function [ weight ] = weight( decode_zero )
%WEIGHT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
l = length(decode_zero);
weight = zeros(1,l/2);
for i=1:l
    if decode_zero(i) == 1
        weight(floor((i+1)/2)) = weight(floor((i+1)/2))+1;
    end
end
for i = 1:l/2
    if weight(i)~=0
        weight(i) = 10*weight(i);
    else
        weight(i)=1;
    end
end

