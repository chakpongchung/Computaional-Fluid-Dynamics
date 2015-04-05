% hwsol1.m -- Used only in hw3plot.m 
% Solves the first order differential equation x'=f(t,x) where the 
% function f(t,x) is as defined by the expression in the string variable 
% ftx. Variables/parameters are obtained from entries in INIT.M, which 
% must be run first. 

global fstring
fstring=ftx;

if t0<ax(2)
 [s,y]=bdode45('eqn1d',t0,ax(2),x0);                  
 plot(s,y,'r') 
 axis(ax) 
 xlabel(['t0 = ',num2str(t0),',   x0 = ',num2str(x0)]) 
end

if ax(1)<t0                                
 [s,y]=bdode45('eqn1d',t0,ax(1),x0);                  
 plot(s,y,'r')
 axis(ax) 
end                                 
                                                 
%hold on                                             
