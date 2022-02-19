clc;clear;close all;
 G=[1 0 0 1 0 1 1 0 ;1 1 1 1 1 0 0 1];
% G = [1 0 0 1 0 0 1 1 0 0;0 1 0 0 1 0 0 1 1 0;1 0 1 0 0 0 0 0 1 1;0 1 0 1 0 1 0 0 0 1];
% G = [1 0 1 1 1 0 1 1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;0 1 0 1 1 0 1 0 0 0 0 0 1 1 1 1;0 1 0 1 0 1 0 1 0 0 1 1 0 0 1 1;0 1 1 0 1 0 0 1 0 1 0 1 0 1 0 1];
t=4;
[k,~]=size(G);
snr_db = 1:0.125:3;
snr = 10.^(0.1*snr_db);
err = 0.5./snr;
input = randi(2,1,1e6)-1;
G_c = [1 1 1;1 0 1];
input_c = code(input,G_c); %ÊäÈëÂë×Ö¾í»ýÂë±àÂë
code = code_q(input_c,G);
ber_cor = zeros(t,length(err));
ber_gkp = zeros(1,length(err));
for m=1:length(err)
code_n = channel(code,err(m));
% code_no = channel(code,0);
% [decode_no,~] = decode_q(code_no,G);
% err_code = find(decode~=decode_no);
[decode,code_zero] = decode_q(code_n,G);
err_zero = [];
for i=1:t
    a=m
    b=i
    temp = err_zero;
    [decode_zero,code_zero] = decode_q(code_zero,G);
    err_zero = find(decode_zero~=0);
    decode_cor = decode;
    for j=1:length(decode_cor)        
            if ~isempty(find(err_zero==j, 1)) || ~isempty(find(temp==j, 1))
                decode_cor(j) = mod(decode_cor(j)+1,2);
            end       
    end
    de_cor = viterbi(G_c,decode_cor);
    ber_cor(i,m) = length(find(de_cor~=input))/length(input);
end
de_gkp = viterbi(G_c,decode);
ber_gkp(m) = length(find(de_gkp~=input))/length(input);
end
ber = [ber_gkp;ber_cor];
%% figure
figure
semilogy(snr_db,ber(1,:))
xlabel('Eb/No')
ylabel('BER')


