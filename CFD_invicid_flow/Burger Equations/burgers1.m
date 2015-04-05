function burgers( method, dx, dt, init_type )
%
%	|-----------------------------------------|
%	| Last Modified Date : 2012.5.22          |
%	| Version : 0.0.1                         |
%	| Author: Yimin Zhong                     |
%	| Platform: Matlab 2009a(Linux)           |
%	|-----------------------------------------|
%
%
%
%Solve Burgers Equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1,
  disp(  'Usage:  burgers( method, dx, dt, init_condition_type )' );
  disp( ['Example: (''dg'', 0.02, 0.015, 3)'] );
  disp( ['Method: ''fdm'' Fridriches']);
  disp(  '        ''fvm'' Modified Fridriches');
  disp(  '        ''dg'' Discontinuous Galerkin with upwind flux');
  disp(  '        ''gv'' Godunov');
  disp(  '        ''upwind'' Upwind Scheme');
  disp( ['Initial Condition: ''1: sine''']);
  disp(  '                 : ''2: triangular''');
  disp(  '                 : ''3: shock''');
  return;
end

if nargin < 2, dx    = 0.02; end;
if nargin < 3, dt    = 0.015; end;
if nargin < 4, init_type= 1;   end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%arguments
tend = 0.24;     % end time
xmax = 1;     % domain is [0,xmax]
nx = floor(xmax/dx);     % spatial discretize
x  = [0 : dx : xmax];    
nt = floor(tend/dt);     % time discretize
dt = tend / nt;
ntprint = 50; % print 50 frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the initial solution values.
u0   = init(x,init_type);
u    = u0;
unew = 0*u;

disp( ['Method: ', method] );
disp( ['   dx = ', num2str(dx)] );
disp( ['   dt = ', num2str(dt)] );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print argument
ntprint = min(nt, ntprint);   
dtprint = tend / ntprint;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step the solution in time.
uall = zeros(ntprint+1,nx+1);
uall(1,:) = u0;
uerr0 = zeros(1,nt);
uerr1 = zeros(1,nt);
uerr2 = zeros(1,nt);

ip = 1;
figure(1)

for i = 1 : nt,

  t = i*dt;
  
  switch lower(method)

    
   case {'fdm','friedriches','a'}
    unew(2:end-1) = 0.5*(u(3:end)+u(1:end-2)) - 0.5*dt/dx * ...
	(f(u(3:end)) - f(u(1:end-2)));
    unew(1)   = 0.5*(u(2)+u(end))-0.5*dt/dx * ...
    (f(u(2))-f(u(end)));
    unew(end) = 0.5*(u(1)+u(end-1))-0.5*dt/dx * ...
    (f(u(1))-f(u(end-1)));
    
   case {'fvm','modified-friedriches','b'}
    unew(2:end-1) = 0.5*(u(3:end)+u(1:end-2))-...
        0.5*(dt/dx)* df(u(2:end-1)).*(u(3:end)-u(1:end-2)); 
	
    unew(1)   = 0.5*(u(2)+u(end))-0.5*dt/dx * u(1).* ...
    (u(2)-u(end));
    unew(end) = 0.5*(u(1)+u(end-1))-0.5*dt/dx * u(end).*...
    (u(1)-u(end-1));
    
   case {'dg','upwind-flux','c'}
    unew(2:end-1) = u(2:end-1) - 0.5*dt/dx * (fn(u(3:end))-fn(u(2:end-1))+...
        fp(u(2:end-1))-fp(u(1:end-2)));
    unew(1) = u(1)- dt/dx * (fn(u(2))-fn(u(1))+fp(u(1))-fp(u(end)));
    unew(end) = u(end) - dt/dx * (fn(u(1))-fn(u(end))+fp(u(end))-fp(u(end-1)));
    
   case {'gv','godunov','d'}
    us = [rpbu(u(1:end-1), u(2:end)), rpbu(u(end),u(1))];
    unew(2:end-1) = u(2:end-1) + dt/dx * (f(us(1:end-2)) - f(us(2:end-1)));
    unew(1)   = u(1)+dt/dx * (f(us(end))-f(us(1)));
    unew(end) = u(end)+dt/dx *(f(us(end-1))-f(us(end)));
   case {'upwind','uw','e'} % upwind (default)
    unew(2:end)= u(2:end) - dt/dx * (f(u(2:end)) - f(u(1:end-1)));
    unew(1)    = u(1)-dt/dx * (f(u(1))-f(u(end)));
    
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot the solution.
  if t >= ip*dtprint,
    plot(x, unew)
    xlabel('x'), ylabel('u')
    title( [method, ' solution at time t=', num2str(t,'%9.4f')] )
    grid on, shg
    pause(0.1)
    ip = ip + 1;
    uall(ip,:) = unew;
  end
  
  uold = u;
  u = unew;

end

figure(2)
nskip = 3;
plot(x,uall(1:nskip:end,:));
xlabel('x'), ylabel('u')
grid on, shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = f( u )
res = 0.5 * u.^2;

function res = df( u )
res = u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = rpbu( uL, uR )
u_star = 0.5 * (uL + uR); 
sig = (uL>=uR);
res =  sig.* (uL .*(u_star>0) + uR .*(u_star<=0)) + (~sig).*(uL .* (uL>0)+ uR .* (uR<0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = fp(x)
res = (max(x,0)).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
function res = fn(x)
res = (min(x,0)).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = init( x, type )

xshift = 0.5;
switch type
 
    case 1     %periodic bc (sine)
        ui = sin(x*pi);
        
    case 2     %periodic bc (triangular)
        ui = 1-2*abs(x-xshift);
        
    case 3     %periodic bc (shock) 
        ui = 1 - 1.5.*(x<=0.25)-1.5*(x>=0.75);
end


