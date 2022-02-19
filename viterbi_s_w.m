function [decoder_output]=viterbi_hard(G,r,weight)
[n,L]=size(G);
n_state=2^((L-1));
for j=0:n_state-1
    for t=0:1
        [state_new,memory_contents]=next_state(j,t,L);
        input(j+1,state_new+1)=t;
        g=rem(memory_contents*G',2);
        newstate(j+1,t+1)=state_new;
         c = char(length(g));
         for i = 1:length(g)
             if g(i) == 1
                  c(i) = '1';
             else
                   c(i) = '0';
             end
         end
        count(j+1,t+1)=bin2dec(c);
    end
end
dist=zeros(n_state,2);
L_r = length(r);
S_length=length(r)/n;
r=reshape(r,n,S_length);
survivor_state=zeros(n_state,S_length+1);
%% 非尾信道解码

for i=1:S_length-L+1
     f=zeros(1,n_state);
    if i<=L
        step=2^(L-i);
    else
        step=1;
    end
    for j=0:step:n_state-1
        for t=0:1
            dist_count=0;
            c_out=double(dec2bin(count(j+1,t+1),n))-48; 
            for m=1:n
                dist_count=dist_count+weight(i)*abs(r(m,i)-c_out(m));
            end
            if ((dist(newstate(j+1,t+1)+1,2)>dist(j+1,1)+dist_count)||f(newstate(j+1,t+1)+1)==0)
                dist(newstate(j+1,t+1)+1,2)=dist(j+1,1)+dist_count;
                survivor_state(newstate(j+1,t+1)+1,i+1)=j;
                f(newstate(j+1,t+1)+1)=1;
            end
        end
    end
    dist=dist(:,2:-1:1);
end
%% 尾信道解码

for i=S_length-L+2:S_length
    f=zeros(1,n_state);
    last_stop=n_state/(2^(i-S_length+L-2));
    for j=0:last_stop-1
        dist_count=0;
        c_out=double(dec2bin(count(j+1,1),n))-48;  
        for m=1:n
            dist_count=dist_count+h_dist(r(m,i),c_out(m));
        end
        if ((dist(newstate(j+1,1)+1,2)>dist(j+1,1)+dist_count)||f(newstate(j+1,1)+1)==0)
            dist(newstate(j+1,1)+1,2)=dist(j+1,1)+dist_count;
            survivor_state(newstate(j+1,1)+1,i+1)=j;
            f(newstate(j+1,1)+1)=1;
        end
    end
    dist=dist(:,2:-1:1);
end

state_sequence=zeros(1,S_length+1);
state_sequence(1,S_length)=survivor_state(1,S_length+1);
for i=1:S_length
   state_sequence(1,S_length-i+1)=survivor_state((state_sequence(1,S_length+2-i)+1),S_length-i+2);
end
decoder_output=zeros(1,S_length-L+1);
for i=1:S_length-L+1
    decoder_output(i)=double(dec2bin(input(state_sequence(1,i)+1,state_sequence(1,i+1)+1),1))-48;
end

