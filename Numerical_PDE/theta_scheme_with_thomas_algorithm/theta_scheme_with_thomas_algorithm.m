% Theta Scheme, with Thomas Algorithm
clear all;
theta =0.6;
alpha = 2; % Control the Standard Deviation
N = 1000; % Steps
T = 2; % End time
t = [0:T/N:T]; 
J = 200; % Space
x = [0:1/J:1];
mu = (T./N)./((1./J).^(2));
U = zeros(N+1,J+1);
mu
theta
% Initial Condition
for j=1:J+1
%============================================================
% Gaussian Curve
 U(1,j) = exp(-alpha.*(x(j)-1./2).^(2));
%============================================================
% Step Function: 
%if (x(j)>0.4)&&(x(j)<0.6)
%    U(1,j) = 1;
%else
%    U(1,j) = 0;
%end
%============================================================
% Boundary Condition
end
for n=1:N+1
    U(n,J+1) = 0;
    U(n,1) = 0;
end
% Implement the Scheme
a = zeros(J+1);
b = zeros(J+1);
c = zeros(J+1);
d = zeros(J+1);
e = zeros(J+1);
f = zeros(J+1);
for n=2:N
    for j=2:J
        a(j) = theta.*mu;
        b(j) = 1+2.*theta.*mu;
        c(j) = theta.*mu;
        d(j) = U(n-1,j)+(1-theta).*mu.*(U(n-1,j+1)-2.*U(n-1,j)+U(n-1,j-1));
        e(j) = c(j)./(b(j)-a(j).*e(j-1));
        f(j) = (d(j)+a(j).*f(j-1))./(b(j)-a(j).*e(j-1)); 
    end
    for k=J:-1:2
        U(n,k) = e(k).*U(n,k+1)+f(k);
    end
end

%Graph plotting

for i=1:20:N+1
    plot(x, U(i,:));
    hold all;
end



