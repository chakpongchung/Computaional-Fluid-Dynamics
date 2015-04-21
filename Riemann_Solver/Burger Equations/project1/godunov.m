clear; 

nx    = 1000; 
 dt    = 0.01; 
 ictype= 5;     % 1 = shock; 2 = expansion; 
                                  % 3 = sonic expansion; 4 = square pulse;5 = sine
                                  
 tend = 2;     % end time
xmax = 16;     % domain length [0,xmax]
dx = xmax/nx; % mesh spacing (constant)
x  = [0 : dx : xmax];
nt = floor(tend/dt);
dt = tend / nt;
ntprint = 50; % for printing              


u0   = uinit(x,ictype);
u    = u0;
unew = 0*u;
us   = unew(1:end-1);

disp( ['   dx = ', num2str(dx)] );
disp( ['   dt = ', num2str(dt)] );

ntprint = min(nt, ntprint);
dtprint = tend / ntprint;

uall = zeros(ntprint+1,nx+1);

uall(1,:) = u0;

ip = 1;
figure(1)
for i = 1 : nt,
  t = i*dt;
  for k=1:size(u,2)-1,
       us(k) = rpbu2(u(k),u(k+1));
  end
 % us = rpbu3(u(1:end-1), u(2:end));
  unew(2:end-1) = u(2:end-1) + dt/dx * (f(us(1:end-1)) - f(us(2:end)));
  unew(1)   = u(1);
  unew(end) = u(end);
    
  % Plot the solution profiles.
  if t >= ip*dtprint,
  %if t <=tend,
    plot(x, unew)
    xlabel('x'), ylabel('u')
    title( ['Solution at time t=', num2str(t,'%9.4f')] )
    grid on, shg
    pause(0.1)
    ip = ip + 1;
    
    uall(ip,:) = unew;
  end
  u = unew;
end


figure(2)
nskip = 3;
plot(x,uall(1:nskip:end,:));
xlabel('x'), ylabel('u')
grid on, shg