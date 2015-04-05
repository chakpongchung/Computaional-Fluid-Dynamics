clear;close all;clc;
N=80;
h=1/N;
X=(1:2*N)'*h/2;
f=X.*(1-X);
e=exp(1);
A=3/(2*e);
B=(4*e-3)/(2*e);
u_exact=A*exp(X)+B*exp(-X)+f-2;
D=zeros(2*N,2*N);
g=zeros(2*N,1);
for i=3:2:2*N-1
    D(i,i)=16/(3*h)+8*h/15;
    D(i,i+1)=-8/(3*h)+h/15;
    D(i+1,i)=-8/(3*h)+h/15;
    D(i,i-1)=-8/(3*h)+h/15;
    D(i-1,i)=-8/(3*h)+h/15;
end
for i=2:2:2*(N-1)
    D(i,i)=14/(3*h)+4*h/15;
    D(i,i+2)=-h/30+1/(3*h);
    D(i+2,i)=-h/30+1/(3*h);
    
end
 D(1,1)=16/(3*h)+8*h/15;
    D(1,2)=-8/(3*h)+h/15;
    D(2,1)=-8/(3*h)+h/15;
    D(2*N,2*N)=7/(3*h)+2*h/15;
    g=zeros(2*N,1);
for j=2:2:2*N
    g(j)=h*f(j)/3;
end
for j=1:2:2*N-1
    g(j)=2*h*f(j)/3;
end
g(2*N)=h/6*f(2*N)-(A*e+B*exp(-1)-2);
C=D\g;
l2norm=norm(C-u_exact,2)
plot([0;X],[0;C],'--O',[0;X],[0;u_exact],'r-');
legend('FEM','u_exact',1);
xlabel('X');
ylabel('u');
title('Quadratic ,N=80');
