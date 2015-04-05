%% Heat equation u_t=u_xx - finite difference scheme - theta method 

%%
%   This program integrates the heat equation  u_t - u_xx = 0
%   on the interval [0,1] using finite difference approximation 
%   via the theta-method. The implicit set of equations are solved at each time step
%   using an LU factorization.


clear all, clc, clf
%% Initial and Boundary conditions
dt = 0.001;
dx = 0.05;

Tf=0.1;
N=floor(Tf/dt);

alpha=20;

% initial conditions
u0 = @(x) 10000000000000000; 

%exp(-alpha.*(x-1./2).^(2));

%1-2*(x-1/2).*sign(x-1/2);

% Dirichlet boundary conditions
left = @(x) 0*x; % left boundary condition
right = @(x) 0*x; % right boundary condition

% the mesh ratio
mu = dt/(dx^2);

tvals=0:dt:Tf;
xvals=0:dx:1;

J=length(xvals);
% Note: the original index j runs from j = 1 ( x = 0) to j = J ( x = 1).

%% Setup of the scheme

% Choose the value of theta:
% theta=0;  % explicit Euler scheme 
% theta=1;  % implicit Euler scheme
% theta=1/2   % Crank-Nicolson scheme

theta =.5;

switch theta
    case 0 
        scheme=' explicit Euler';
    case 1
        scheme=' implicit Euler';
    case 1/2
        scheme=' Crank-Nicolson';
end

u=zeros(J,N);

u(:,1)=u0(xvals);

E = ones(J ,1);
D = spdiags([-E 2*E -E],[-1,0,1],J,J);
I = speye(J);

A = I+theta*mu*D;
B = I-(1- theta )*mu*D;

A(1,:) = 0; A(1,1) = 1;
A(J,:) = 0; A(J,J) = 1;

[L,U] = lu(A);

%% Time iteration

Nframes=5;
heatmovie=moviein(Nframes);

snap0=u0(xvals);

n=0;
nframe=floor(N/Nframes);
N=nframe*Nframes;

for k=1:Nframes;
for m = 1:nframe
    n=n+1; % counter for iteration
    rhs=B*u(:,n);
    rhs([1,J])=0;
    y = L\rhs;
    u(:,n+1) = U\y;
end 

snap1 = [left(n*dt);u(2:J-1,n);right(n*dt)];
lg=sprintf(' time t=%g \n',n*dt);
pause(0.5)

plot(xvals,snap0,'--b',xvals,snap1, '--mo')
xlabel('x'), ylabel('u')
legend('initial condition','computed solution')
title(strcat('Numerical solution at ', lg,'using ',scheme,', \mu= ',num2str(mu)));
hold on

heatmovie(:,k)=getframe;
end
title(strcat('Numerical solution at different times using ',scheme,', \mu= ',num2str(mu)));

%% Plot the final results

hold off
plot(xvals,snap0,'--b',xvals, snap1,'-ro'),
legend('initial condition','computed solution')

title(strcat('Numerical solution at ', lg,' using ',scheme,', \mu= ',num2str(mu)));


figure, surf(xvals,tvals,u'), xlabel x; ylabel t; zlabel u,  rotate3d on
