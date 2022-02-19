function [ code_out ] = code_q( input_c,G )
%CODE_Q 此处显示有关此函数的摘要
%   此处显示详细说明
[c,n] = size(G);
n = n/2;
k = n-c;
l = length(input_c);
input_cin = reshape(input_c',[k,l/k])';%量子卷积码电路每次输入k个
r = rank(G(1:c,1:c));
A1 = G(1:r,r+1:n-k);
C1 = G(1:r,n+r+1:2*n-k);
A2 = G(1:r,n-k+1:n);
C2 = G(1:r,2*n-k+1:2*n);
B = G(1:r,n+1:n+r);
D = G(r+1:n-k,n+1:n+r);
E = G(r+1:n-k,2*n-k+1:2*n);
X = mod([zeros(k,r),E',eye(k),E'*C1'+C2',zeros(k,n-r)],2);
Z = [A2',zeros(k,n-k-r),eye(k)];
T = [];
u2 = X(1:k,r+1:n-k);
s = [];
code_out = [];
if r<n-k
    for m = 1:k
        temp = zeros(1,n);
        st = [];
        for i = r+1:n-k
            temp(i) = 1;
            if i == r+1               
                st = u2(m,i-r)*temp;
            else
                st = st+u2(m,i-r)*temp;
            end
        end
        s = [s;st];
    end
end
for ii = 1:l/k
        code_in = [zeros(1,n-k),input_cin(ii,:)];
        if r<n-k
            for i = 1:k
                flag = find(s(i,:)==1);
                for j = 1:length(flag)
                    if code_in(n-k+i)==1
                        code_in(flag(j)) = mod(code_in(flag(j))+1,2);
                    end
                end
            end
        end
        for i=1:n-k
            temp = G(i,:);
            temp(i) = 0;
            temp(i+n) = 0;
            T = [T;temp];
        end
        code = code_in;
        for i=1:c
            code = [code;code];
            code(2^(i-1)+1:2^i,i) = 1;
            for j = 2^(i-1)+1:2^i
                for m = n+1:2*n
                    if T(i,m) == 1 && code(j,m-n) == 1
                        code(j,m-n) = -code(j,m-n);
                    end
                end
                for kk = 1:n
                    if T(i,kk) == 1
                        code(j,kk) = ~code(j,kk);
                    end
                end
            end
        end
code_out = [code_out;code];
end
end

