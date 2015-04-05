function ploterr
h=zeros(1,4);
l2norm=h;
N=[10 20 40 80];
for i=1:4
    [h(i),l2norm(i)]=linearElement(N(i));
end
loglog(h,l2norm,'o')
end
function [h,l2norm]=linearElement(N)
%clear;close all;clc;
%N=10;

if nargin<1, N = 10; end
h=1/N;
X=(1:N)'*h;
f=X.*(1-X);
e=exp(1);
A=3/(2*e);
B=(4*e-3)/(2*e);
u_exact=A*exp(X)+B*exp(-X)+f-2;
D=zeros(N,N);
for i=1:N-1
    D(i,i)=2/h+2*h/3;
    D(i,i+1)=-1/h+h/6;
    D(i+1,i)=-1/h+h/6;
end
D(N,N)=1/h+h/3;
g=h*f;
g(N)=-(A*e+B*exp(-1)-2);
C=D\g;
l2norm=norm(C-u_exact,2);
end



