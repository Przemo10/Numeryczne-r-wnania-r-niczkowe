
function [s,U,t]=linshoot17(p,q,f,a,b,al,be,N)
%solves BVP -u''+p(t)u'+q(t)u=f(t) with u(a)=al u(b)=be
%using shooting method i.e. solves IVP
%-u''+p(t)u'+q(t)u=f(t) with u(a)=al u'(a)=s for s=0,1
%getting u(T;0) and u(T;1) then u(t;s) with
% s=(be - u(T;0))/(u(T;1)-u(T;0) 
%solves BVP 
%INPUT 
%p,q,f - functio handles to parameter functions
%a,b - boundary points
%al,be bnd values
%N - no of discrete time points i.e. t(k)=a+k*h h=(b-a)/N k=0,..,N
%OUTPUT
%s - shooting value i.e .u'(a)
%U N+1 x 2 matrix with discrete solutions i.e. U(:,k) approx of [u(t(k),u'(t(k)]
%EXAMPLE: p=q=f=@(t) 0.0;q=@(t) -1.0; a=0;b=1;al=be=1;N=100;[s,U,t]=linshoot17(p,q,f,a,b,al,be,N);U(end,1)-be
%same but b=10;%b=10;[s,U,t]=linshoot17(p,q,f,a,b,al,be,N);U(end,1)-be
global pp qq ff
pp=p;qq=q;ff=f;
t=linspace(a,b,N+1);

U1=lsode(@vf,[al;0],[a,b]);
U2=lsode(@vf,[al;1],[a,b]);
ga1=U1(2,1);
ga2=U2(2,1);
if(abs(ga2-ga1)<1e-14)
 printf("Warning: bnd may have no solutions!\n");
end
s=(be-ga1)/(ga2-ga1);
U=lsode(@vf,[al;s],t);

end

function y=vf(x,t)
global ff pp qq
A=[0,1;-qq(t),-pp(t)];
y=A*x+[0;ff(t)];
end
