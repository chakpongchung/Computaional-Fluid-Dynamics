%   CHUNG,Chak Pong


clear all, clc, clf;

alpha=0.05;
J=20;
x=[0:1/J:1];

T=100;
t=[0:1/T:1];

U=zeros(T+1,J+1);

% Dirichlet boundary conditions, ininitialed as 0;



%I.C. initialization
for j=2:J

 U(1,j) = exp(-alpha.*(x(j)-0.5).^(.5));
 
end

% the mesh ratio as shown in page 24;also the Mu and theta
Mu=(1/T)/((1/J)^2);


%% Setup of the scheme
theta=0.7;

a=zeros(1);
b=zeros(1);
c=zeros(1);
d=zeros(1);
e=zeros(1);
f=zeros(1);

%notation of the coefficient follows the book in page 24,equation 2.65

for j=1:J+1
    a(j)=theta*Mu;
    b(j)=1+2*theta*Mu;
    c(j)=theta*Mu;
    
    
end;

for j=2:J+1
    e(j)=c(j)/(b(j)-a(j)*e(j-1));
end
 
for n=1:T
    for j=2:J
        d(n,j)=U(n,j)+[(1-theta)*Mu]*[U(n,j+1)-2*U(n,j)+U(n,j-1)];
        f(j)=[d(n,j)+a(j)*f(j-1)]/[b(j)-a(j)*e(j-1)];
    end
    
 %Now using Thomas Algorithm to do backward substitution   
    for j=J:-1:2
        U(n+1,j)=e(j)*U(n+1,j+1)+f(j);
    end
 %backward substitution done for time step n 
end


%plot3

%% plot

for i=1:T+1
    plot(x, U(i,:));
    hold all;
end







