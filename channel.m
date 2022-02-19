function [ out ] = channel( code,snr )
%CHANNEL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[c,r] = size(code);
out = zeros(c,r);
noise = [];
for i=1:c
    for j=1:r
        err =  snr*0.625*randn(1);
        cor = 0.3*0.625*snr*randn(1);
        cor = mod(cor+err,1);
        if cor>0.5
            cor = cor-1;
        end
        cor = err-cor;
        out(i,j) = mod(1+code(i,j) + cor,2)-1 ;
        noise = [noise,abs(err)];
    end
end
noise = mean(noise)
end

