clc;clear;
%generator matrix
G=[1 0 0 1 0 1 1 0 ;1 1 1 1 1 0 0 1];
input = [1 1];
%% code
G = abs(rref(G));
[c,r] = size(G);
n = 4;
k = 2;
l = length(input);
input = reshape(input',[k,l/k])';
r = rank(G);
A1 = G(1:r,r+1:n-k);
C1 = G(1:r,n+r+1:2*n-k);
A2 = G(1:r,n-k+1:n);
C2 = G(1:r,2*n-k+1:2*n);
B = G(1:r,n+1:n+r);
D = G(r+1:n-k,n+1:n+r);
E = G(r+1:n-k,2*n-k+1:2*n);
X = [zeros(k,r),E',eye(k),E'*C1'+C2',zeros(k,n-r)];
Z = [A2',zeros(k,n-k-r),eye(k)];
T = [];
decode_final = [];
for ii = 1:l/k
    code_in = [zeros(1,n-k),input(ii,:)];
    for i=1:c
        temp = G(i,:);
        temp(i) = 0;
        temp(i+n) = 0;
        T = [T;temp];
    end
    code = code_in;
    sign = 1;
    for i=1:c
        code = [code;code];
        sign = [sign,sign];
        code(2^(i-1)+1:2^i,i) = 1;
        for j = 2^(i-1)+1:2^i
            for m = n+1:2*n
                if T(i,m) == 1 && code(j,m-n) == 1
                    sign(j) = -sign(j);
                end
            end
            for kk = 1:n
                if T(i,kk) == 1
                    code(j,kk) = ~code(j,kk);
                end
            end
        end
    end
    %% channel\
    code(:,4) = [1 1 1 1]';
    
    
    
    %% decode
    decode = [code';zeros(k,2^c)]';
    for t = 1:k
        for i=1:2^c
            for j=1:n
                if Z(t,j)==1 && decode(i,j)==1
                    decode(i,n+t) = ~decode(i,n+t);
                end
            end
            for j=n+1:2*n
                if X(t,j)==1 && decode(i,n+t)==1 && decode(i,j-n)==1
                    sign(i) = -sign(i);
                end
            end
            for j=1:n
                if X(t,j)==1 && decode(i,n+t)==1 
                    decode(i,j) = ~decode(i,j);
                end
            end
        end
    end
    decode_last = decode(1,n+1:end);
    decode_final = [decode_final,decode_last]
end