% burgers.m: Solve the Burgers' equation using a variety of methods. 
%
% Date:   March 5, 2011
% Author: John Stockie
%         Department of Mathematics
%         Simon Fraser University
%
% Method choices: 'Upwind/uw', 'Upwind(non-cons.)/uwnc', 
%                 'Lax-Friedrichs/lf', 'Lax-Wendroff/lw', 
%                 'MacCormack/mc'.
%
% Initial conditions: 1 = piecewise constant (shock)
%                     2 = piecewise constant (expansion)
%                     3 = cosine
%
% Examples:
% 1. Try Lax-Wendroff (ictype=1) with nx=100, dt=0.0105->0.011

function burgers( method, nx, dt, ictype )

if nargin < 1,
  disp(  'Usage:  burgers( method, nx, dt, ictype )' );
  disp( ['        where method is one of ''Upwind/uw'', ']);
  disp( ['        ''Upwind(non-cons.)/uwnc'', ''Lax-Friedrichs/lf'',']);
  disp( ['        ''Lax-Wendroff/lw'', ''MacCormack/mc'''] );
  disp( ['Defaults: (''upwind'', 100, 0.1, 1)'] );
  return;
end

if nargin < 2, nx    = 100; end;
if nargin < 3, dt    = 0.1; end;
if nargin < 4, ictype= 1;   end;  % 1 = shock; 2 = expansion; 3 = cos

tend = 4;     % end time
xmax = 8;     % domain is [0,xmax]

dx = xmax/nx % mesh spacing (constant)
x  = [0 : dx : xmax];
nt = floor(tend/dt);
dt = tend / nt;
ntprint = 50; % for printing

% Set up the initial solution values.
u0   = uinit(x,ictype);
u    = u0;
unew = 0*u;

disp( ['Method: ', method] );
disp( ['   dx = ', num2str(dx)] );
disp( ['   dt = ', num2str(dt)] );

ntprint = min(nt, ntprint);
dtprint = tend / ntprint;

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
    
   case {'upwind(non-cons.)','uwnc'}
    unew(2:end) = u(2:end) - dt/dx * u(2:end) .* (u(2:end) - u(1:end-1));
    unew(1)   = u(1);
    
   case {'lax-friedrichs','lf'}
    unew(2:end-1) = 0.5*(u(3:end)+u(1:end-2)) - 0.5*dt/dx * ...
	(f(u(3:end)) - f(u(1:end-2)));
    unew(1)   = u(1);
    unew(end) = u(end);
    
   case {'lax-wendroff','lw'}
    unew(2:end-1) = u(2:end-1) ...
	- 0.5*dt/dx * (f(u(3:end)) - f(u(1:end-2))) ...
	+ 0.5*(dt/dx)^2 * ... 
	( df(0.5*(u(3:end) + u(2:end-1))) .* (f(u(3:end)) - f(u(2:end-1))) - ...
	  df(0.5*(u(2:end-1) + u(1:end-2))) .* (f(u(2:end-1)) - f(u(1:end-2))) );
    unew(1)   = u(1);
    unew(end) = u(end);
    
   case {'maccormack','mc'}
    us           = u(1:end-1) - dt/dx * (f(u(2:end)) - f(u(1:end-1)));
    unew(2:end-1)= 0.5*(u(2:end-1) + us(2:end)) - ...
	0.5*dt/dx* (f(us(2:end)) - f(us(1:end-1)));    
    unew(1)   = u(1);
    unew(end) = u(end);
    
   otherwise  % upwind (default)
    unew(2:end)= u(2:end) - dt/dx * (f(u(2:end)) - f(u(1:end-1)));
    unew(1)    = u(1);
    
  end

  % Plot the solution profiles.
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

% Plot a sequence of profiles throughout the run.
figure(2)
nskip = 3;
plot(x,uall(1:nskip:end,:));
xlabel('x'), ylabel('u')
grid on, shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret = f( u )
ret = 0.5 * u.^2;

function ret = df( u )
ret = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = uinit( x, ictype )

xshift = 1.0;
switch ictype
 case 1    % shock (uL > uR)
  uL = 1.0;
  uR = 0.5;
  ui = uR + (uL-uR) * ((x-xshift) <= 0.0);
 case 2    % expansion (uL < uR)
  uL = 0.5;
  uR = 1.0;
  ui = uR + (uL-uR) * ((x-xshift) <= 0.0);
  
      case 3
        ui = ( sin(0.25*pi*(x-0)) );
 otherwise % cosine (forms shock later)
  ui = 0.5*(1 + cos(pi*(x-xshift))) .* (abs(x-xshift) <= 1);
end


