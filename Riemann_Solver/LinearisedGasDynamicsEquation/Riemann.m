% Riemann solver for linear dynamic equation
% ul - 1x2 array, state on the left
% ur - 1x2 array, state on the right
% flux - 1x2 array, output flux
% Author Luo Li, lluoac@ust.hk
function [flux]=Riemann(ul,ur)
% parameters
a=1.0;
rho0=1.0;
flux=zeros(1,2);
alpha2=(a*ul(1)+rho0*ul(2))/(2.0*a*rho0);
beta1=(a*ur(1)-rho0*ur(2))/(2.0*a*rho0);
u1=beta1*rho0+alpha2*rho0;
u2=beta1*(-a)+alpha2*a;
flux(1,1)=rho0*u2;
flux(1,2)=a*a/rho0*u1;


