clear;clc;close all;

num = 50;
t=num;

% number of space gridpoints
dt = 1/num;                                          % time step
t_final = 1;                                            % final time
a = 1;                                  
Q=zeros(num+1,num+1);

dx = 1/(num);
x=0:dx:1;

mu=a*dt/dx;



for i=1:num


Q(1,i)=exp(-0.1*(i-num/20).^2);
%square initial condition
%if (x(i)>0.4)&&(x(i)<0.6)
%    Q(1,i) = 11;
%else
%    Q(1,i) = 0;
end
%============================================================
% Boundary Condition



for n=1:num
    Q(n,1)= 0;
    Q(n,num)= 0;
end

limiter=Q;        %%here comes to the limiter !!!!

for t=1:num
    for i=2:num
        
        limiter(t,i)= (  Q(t,i)-Q(t,i-1)  )/ (2 * dx);  %% As you may notice,it has some problem on for i =1,so I just start from i =2
        
                                                         %%This is Beam-Warming limiter,you can check the book,page 107,for minmod,check page 111
    end
end


for t=1:num
   for i=2:num
      %%Q(t+1,i)=mu*Q(t,i-1) + (1-mu)*Q(t,i);
      
      Q(t+1,i)=Q(t,i)-mu*(  Q(t,i)-Q(t,i-1) )- 0.5*mu *(dx -a*dt) * (limiter(t,i)-limiter(t,i-1) ) 
   end
   
    plot(x,Q(t,:),'--rs','LineWidth',2,...    %%at the 20th time level,you can change this from 1 to 100
                 'MarkerEdgeColor','k',...
                 'MarkerFaceColor','g',...
                 'MarkerSize',2);
    axis([0 1.5 -0.1 1])
    title(sprintf('t=%0.2f',t*dt))
    
    %%pause(0.0001) ;
    drawnow
end

%%hold on;
space_step=1:1:100;

time_interval=1:1:100;
    
% plot(space_step,Q(20,time_interval),'--rs','LineWidth',2,...    %%at the 20th time level,you can change this from 1 to 100
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',2);

                
            






