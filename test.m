clc;clear;close all;
%generator matrix
G=[1 0 0 1 0 1 1 0 ;1 1 1 1 1 0 0 1];
% G = [1 0 0 1 0 0 1 1 0 0;0 1 0 0 1 0 0 1 1 0;1 0 1 0 0 0 0 0 1 1;0 1 0 1 0 1 0 0 0 1];
% G = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;0 1 0 1 1 0 1 0 0 0 0 0 1 1 1 1;0 1 0 1 0 1 0 1 0 0 1 1 0 0 1 1;0 1 1 0 1 0 0 1 0 1 0 1 0 1 0 1];
ber_gkp = [];
ber_gkp_s = [];
ber_gkp_w = [];
err = 0.05:0.0125:0.3;
input = randi(2,1,1000)-1;
%% code 量子卷积码编码
G_c = [1 1 1;1 0 1];  %传统卷积码的生成矩阵
input_c = code(input,G_c); %输入码字卷积码编码
G = mod(abs(rref(G)),2);%量子卷积码生成矩阵的标准形式
n = 4;k = 2;
c = n-k;
%n = 8;k = 3;
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
for num = 1:length(err)
    weight_final = [];
    decode_final_gkp = [];
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

       %% channel 信号过高斯位移信道
        code_gkp = code;
        gkp_err = zeros(2.^c,n);
        gkp_cor = zeros(2.^c,n);
        for i = 1:2.^c
            for j = 1:n                            
                gkp_err(i,j) =  err(num)*randn(1);
                gkp_cor(i,j) = 0.3*err(num)*randn(1); 
            end
        end
        %% gkp decode  
        [g_c,g_r] = size(gkp_cor);
        for ii = 1:g_c
            for jj = 1:g_r
                gkp_cor(ii,jj) = mod(gkp_cor(ii,jj)+gkp_err(ii,jj),1);
                if gkp_cor(ii,jj)>0.5
                    gkp_cor(ii,jj) = gkp_cor(ii,jj)-1;
                end
                gkp_cor(ii,jj) = gkp_err(ii,jj)-gkp_cor(ii,jj);
            end
        end
        gkp_out = mod(1+code_gkp + gkp_cor,2)-1 ;
        decode_gkp = [gkp_out';zeros(k,2^c)]';
        for t = 1:k
            for i=1:2^c
                for j=1:n
                    if Z(t,j)==1 
                        if decode_gkp(i,n+t)>0
                            decode_gkp(i,n+t) = decode_gkp(i,j).^2*(1-decode_gkp(i,n+t))+sqrt(1-decode_gkp(i,j).^2)*decode_gkp(i,n+t);
                        else
                            decode_gkp(i,n+t) = -decode_gkp(i,j).^2*(1+decode_gkp(i,n+t))+sqrt(1-decode_gkp(i,j).^2)*decode_gkp(i,n+t);
                        end
                    end
                end
            end
        end
        for  t=1:k
            for i = 1:2^c
                for j=n+1:2*n
                    if X(t,j)==1 
                        decode_gkp(i,j-n) = -decode_gkp(i,n+t).^2*decode_gkp(i,j-n) + (1-decode_gkp(i,n+t).^2)*decode_gkp(i,j-n);
                    end
                end
                for j=1:n
                    if X(t,j)==1 
                        if decode_gkp(i,j)>0
                            decode_gkp(i,j) = decode_gkp(i,n+t).^2*(1-decode_gkp(i,j))+sqrt(1-decode_gkp(i,n+t).^2)*decode_gkp(i,j);
                        else
                            decode_gkp(i,j) = -decode_gkp(i,n+t).^2*(1+decode_gkp(i,j))+sqrt(1-decode_gkp(i,n+t).^2)*decode_gkp(i,j);
                        end
                    end
                end
            end
        end
        decode_last_gkp = decode_gkp(1,n+1:end);
        decode_final_gkp = [decode_final_gkp,decode_last_gkp];
 %% error weight gkp
        decode_zero = abs(decode_gkp(1:n,1:n));
        for i=k:-1:1
            for j=n:-1:k+1
                for nn=1:n
                    if   T(i,j)==1
                        decode_zero(nn,j)=decode_zero(nn,i)*abs(abs(decode_zero(nn,j))-1);
                    end
                end
            end
           for nn=1:n
               decode_zero(nn,i)=0;
           end
        end
        for i=1:n
            for j=1:n
                if decode_zero(i,j)<0.5
                    decode_zero(i,j)=0;
                else
                    decode_zero(i,j)=1;
                end
            end
        end
        weight = length(find(decode_zero)~=0);
        if weight == 0
            weight = 1;
        else
            weight = 100*weight;
        end
        weight_final = [weight_final,weight];
    end
    decode_final_gkp = abs(decode_final_gkp);
    %% gkp conv
    decode_final_s = decode_final_gkp;
    for i = 1:length(decode_final_gkp)
        if decode_final_gkp(i)>0.5
            decode_final_gkp(i) = 1;
        else
            decode_final_gkp(i) = 0;
        end
    end
    decode_final_gkp_w = decode_final_gkp;
    for i = 1:length(decode_final_gkp)
        if weight_final(floor((i+1)/2))>1
            decode_final_gkp_w(i) = mod(decode_final_gkp_w(i)+1,2);
        end
    end
    de_gkp = [];
    de_gkp = viterbi(G_c,decode_final_gkp);
    ber_gkp = [ber_gkp,length(find(de_gkp~=input))/length(input)];
    de_gkp_w = [];
    de_gkp_w = viterbi_s_w(G_c,decode_final_gkp,weight_final);
    ber_gkp_w = [ber_gkp_w,length(find(de_gkp_w~=input))/length(input)];
    de_gkp_s = [];
    de_gkp_s = viterbi_s(G_c,decode_final_s);
    ber_gkp_s = [ber_gkp_s,length(find(de_gkp_s~=input))/length(input)];
end
%% test
weight_diff = find(weight_final~=1);
diff = find(input_c~=decode_final_gkp);
%% figure
figure
semilogy(err,ber_gkp,'b-*',err,ber_gkp_w,'r',err,ber_gkp_s,'p')
xlabel('Eb/No')
ylabel('BER')
legend('o','w','s')

