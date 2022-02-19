clear;clc;close all;
input = randi(2,1,1e6)-1;
G_c = [1 1 1;1 0 1]; 
c = code(input,G_c);
c = awgn(c,10);
for i = 1:length(c)
    if c(i)>0.5
        c(i)=1;
    else
        c(i)=0;
    end
end
decode = Copy_of_viterbi(G_c,c);
ber = length(find(decode~=input))/length(input);
