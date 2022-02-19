%£¨n£¬1£¬L£©¾í»ıÂë±àÂë
 function [output]=code(input,g)
n=length(input);
li=size(g,2); 
n0=size(g,1); 
u=[zeros(size(1:li-1)),input,zeros(size(1:(li-1)))];
u1=u(li:-1:1);
for i=1:n+li-2
   u1=[u1,u((i+li):-1:i+1)];
end
uu=reshape(u1,li,n+li-1);
G=rem(g*uu,2);
output=reshape(rem(g*uu,2),1,n0*(n+li-1));
 end