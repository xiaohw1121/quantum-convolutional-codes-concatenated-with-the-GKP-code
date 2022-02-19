function [ decode_final_gkp,code_zero] = decode_q(code,G )
%DECODE_Q 此处显示有关此函数的摘要
%   此处显示详细说明
[c,n] = size(G);
n = n/2;
k = n-c;
r = rank(G(1:c,1:c));
C1 = G(1:r,n+r+1:2*n-k);
A2 = G(1:r,n-k+1:n);
C2 = G(1:r,2*n-k+1:2*n);
E = G(r+1:n-k,2*n-k+1:2*n);
X = mod([zeros(k,r),E',eye(k),E'*C1'+C2',zeros(k,n-r)],2);
Z = [A2',zeros(k,n-k-r),eye(k)];
[l,~] = size(code);
l = l/2^c;
decode_final_gkp = [];
code_zero = [];
for num = 1:l
    gkp_out = code(2^c*num-2^c+1:2^c*num,1:n);
    decode = [gkp_out';zeros(k,2^c)]';
    for t = 1:k
        for i=1:2^c
            for j=1:n
                if Z(t,j)==1
                    if decode(i,n+t)>0
                        decode(i,n+t) = decode(i,j).^2*(1-decode(i,n+t))+sqrt(1-decode(i,j).^2)*decode(i,n+t);
                    else
                        decode(i,n+t) = -decode(i,j).^2*(1+decode(i,n+t))+sqrt(1-decode(i,j).^2)*decode(i,n+t);
                    end
                end
            end
        end
    end
    for  t=1:k
        for i = 1:2^c
            for j=n+1:2*n
                if X(t,j)==1
                    decode(i,j-n) = -decode(i,n+t).^2*decode(i,j-n) + (1-decode(i,n+t).^2)*decode(i,j-n);
                end
            end
            for j=1:n
                if X(t,j)==1
                    if decode(i,j)>0
                        decode(i,j) = decode(i,n+t).^2*(1-decode(i,j))+sqrt(1-decode(i,n+t).^2)*decode(i,j);
                    else
                        decode(i,j) = -decode(i,n+t).^2*(1+decode(i,j))+sqrt(1-decode(i,n+t).^2)*decode(i,j);
                    end
                end
            end
        end
    end
    decode_last_gkp = decode(1,n+1:end);
    code_zero = [code_zero;decode(1:2^c,1:n)];
    decode_final_gkp = [decode_final_gkp,decode_last_gkp];
end
decode_final_gkp = abs(decode_final_gkp);
for i = 1:length(decode_final_gkp)
    if decode_final_gkp(i)>0.5
        decode_final_gkp(i) = 1;
    else
        decode_final_gkp(i) = 0;
    end
end
end

