% Solve linear dynamic equation by exact Riemann solver and SSP-RK method
% n - number of spatial grid in one-dimension
% cfl - CFL number
% t_F - final time
% bc -  boundary condition type: 1 - open boundary, 2 - reflection boundary
% flg - indicate number for different cases: 1 - square perturbation wave, 2 - triangle
% perturbation wave, 3 - sin wave
% called at MATLAB command line:  linear_dynamic(500, 0.7, 1.0, 1, 2)
% Author: Luo Li, lluoac@ust.hk
function linear_dynamic(n,cfl,t_F,bc,flg)
itmax=1000; % maximum of time steps
m=2; % number of ghost points
xmin=-0.5;
xmax=0.5;
dx=(xmax-xmin)/(n-1); % space step length
x=zeros(n,1);
u0=zeros(n+2*m,2);  % RK 0
u1=zeros(n+2*m,2);  % RK 1
u2=zeros(n+2*m,2);  % RK 2
u3=zeros(n+2*m,2);  % RK 3
% parameters
a=1.0;
% initialize
if flg==1 % square perturbation wave
    for i=1:n
        x(i)=-0.5+(i-1)*dx;
        if x(i)<-0.1
           u0(m+i,1)=0.0;
           u0(m+i,2)=0.0;
        elseif x(i)<=0.1
           u0(m+i,1)=0.1;
           u0(m+i,2)=0.0;
        else
           u0(m+i,1)=0.0;
           u0(m+i,2)=0.0;
        end
    end
elseif flg==2 % triangle perturbation wave
    for i=1:n
        x(i)=-0.5+(i-1)*dx;
        if x(i)<-0.05
           u0(m+i,1)=0.0;
           u0(m+i,2)=0.0;
        elseif x(i)<=0.0
           u0(m+i,1)=sqrt(3.0)/2.0*(x(i)+0.05);
           u0(m+i,2)=0.0;
        elseif x(i)<=0.05
           u0(m+i,1)=-sqrt(3.0)/2.0*(x(i)-0.05);
           u0(m+i,2)=0.0;            
        else
           u0(m+i,1)=0.0;
           u0(m+i,2)=0.0;
        end
    end 
elseif flg==3 % sin wave
    for i=1:n
        x(i)=-0.5+(i-1)*dx;
        if x(i)<-0.1
           u0(m+i,1)=0.0;
           u0(m+i,2)=0.0;
        elseif x(i)<=0.1
           u0(m+i,1)=0.1*sin(pi*(10.0*x(i)+1.0));
           u0(m+i,2)=0.0;
        else
           u0(m+i,1)=0.0;
           u0(m+i,2)=0.0;
        end
    end
end
if bc==1
% open boundary
for i=1:m
    u0(i,:)=u0(m+1,:);
    u0(m+n+i,:)=u0(m+n,:);
end
elseif bc==2 
% reflection boundary
for i=1:m
    % left 
    u0(i,1)=u0(2*m+1-i,1);  % density not change
    u0(i,2)=-u0(2*m+1-i,2); % velocity changes
    % right
    u0(m+n+i,1)=u0(m+n+1-i,1);
    u0(m+n+i,2)=-u0(m+n+1-i,2);
end    
end
subplot(2,1,1)
plot(x,1.0+u0(m+1:m+n,1),'LineWidth',2);
axis([-0.5 0.5 0.7 1.3]); % called after plot()
subplot(2,1,2)
plot(x,u0(m+1:m+n,2),'--','LineWidth',2);
axis([-0.5 0.5 -0.3 0.3]); % called after plot()

% start time stepping
t=0.0;
for it=1:itmax
    % compute time step length dt=cfl*dx/|a|
    dt=cfl*dx/abs(a);
    % final time step length
    if t<t_F && t+dt>t_F
        dt=t_F-t;
    end
    t=t+dt
    pause;
    % if time is up, stop time stepping
    if t+dt>t_F || it==itmax
        break;
        subplot(2,1,1)
        plot(x,1.0+u0(m+1:m+n,1),'LineWidth',2);
        axis([-0.5 0.5 0.7 1.3]); % called after plot()
        subplot(2,1,2)
        plot(x,u0(m+1:m+n,2),'--','LineWidth',2);
        axis([-0.5 0.5 -0.3 0.3]); % called after plot()
    end
    
    lambda=dt/dx;
    % RK step 1
    for i=1:n
        u1(m+i,:)=u0(m+i,:)-lambda*(Riemann(u0(m+i,:),u0(m+i+1,:))-Riemann(u0(m+i-1,:),u0(m+i,:)));
    end
    if bc==1
    % open boundary
    for i=1:m
        u1(i,:)=u1(m+1,:);
        u1(m+n+i,:)=u1(m+n,:);
    end
    elseif bc==2 
    % reflection boundary
    for i=1:m
        u1(i,1)=u1(2*m+1-i,1);    
        u1(i,2)=-u1(2*m+1-i,2);
        u1(m+n+i,1)=u1(m+n+1-i,1);
        u1(m+n+i,2)=-u1(m+n+1-i,2);
    end    
    end
    % RK step 2
    for i=1:n
        u2(m+i,:)=3.0/4.0*u0(m+i,:)+1.0/4.0*(u1(m+i,:)-lambda*(Riemann(u1(m+i,:),u1(m+i+1,:))-Riemann(u1(m+i-1,:),u1(m+i,:))));
    end
    if bc==1
    % open boundary
    for i=1:m
        u2(i,:)=u2(m+1,:);
        u2(m+n+i,:)=u2(m+n,:);
    end
    elseif bc==2 
    % reflection boundary
    for i=1:m
        u2(i,1)=u2(2*m+1-i,1);    
        u2(i,2)=-u2(2*m+1-i,2);
        u2(m+n+i,1)=u2(m+n+1-i,1);
        u2(m+n+i,2)=-u2(m+n+1-i,2);
    end    
    end 
    % RK step 3
    for i=1:n
        u3(m+i,:)=1.0/3.0*u0(m+i,:)+2.0/3.0*(u2(m+i,:)-lambda*(Riemann(u2(m+i,:),u2(m+i+1,:))-Riemann(u2(m+i-1,:),u2(m+i,:))));
    end
    if bc==1
    % open boundary
    for i=1:m
        u3(i,:)=u3(m+1,:);
        u3(m+n+i,:)=u3(m+n,:);
    end
    elseif bc==2 
    % reflection boundary
    for i=1:m
        u3(i,1)=u3(2*m+1-i,1);    
        u3(i,2)=-u3(2*m+1-i,2);
        u3(m+n+i,1)=u3(m+n+1-i,1);
        u3(m+n+i,2)=-u3(m+n+1-i,2);
    end    
    end 
    
    % swap
    u0(:,:)=u3(:,:);
    if bc==1
    % open boundary
    for i=1:m
        u0(i,:)=u0(m+1,:);
        u0(m+n+i,:)=u0(m+n,:);
    end
    elseif bc==2 
    % reflection boundary
    for i=1:m
        u0(i,1)=u0(2*m+1-i,1);    
        u0(i,2)=-u0(2*m+1-i,2);
        u0(m+n+i,1)=u0(m+n+1-i,1);
        u0(m+n+i,2)=-u0(m+n+1-i,2);
    end    
    end
    subplot(2,1,1)
    plot(x,1.0+u0(m+1:m+n,1),'LineWidth',2);
    axis([-0.5 0.5 0.7 1.3]); % called after plot()
    subplot(2,1,2)
    plot(x,u0(m+1:m+n,2),'--','LineWidth',2);
    axis([-0.5 0.5 -0.3 0.3]); % called after plot()
end
