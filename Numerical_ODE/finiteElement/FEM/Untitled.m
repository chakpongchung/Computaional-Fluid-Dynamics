% Implementation of Midpoint, Trapezoidal and Simpson Rule.
% Mohammad Ismail Hossain 
% jacobs University Bremen

clc
clear all;

a = 5;
b = 10;
c = 4;

Int_mdp = zeros(1,100);
Int_trpz = zeros(1,100);
Int_simp = zeros(1,100);
for N = 1:100
    
    % Midpoint Rule
    dx = (b-a)/N;
    sum_mdp = 0;
    for p = 1:N
        X = a + p*(dx)-(dx/2); 
        Y = sqrt(((X^2)-(c^2)))/X;
        sum_mdp = sum_mdp + Y*dx;      
    end
    Int_mdp(N) = sum_mdp;
    
    % Trapezoidal Rule
    dx = (b-a)/(N-1);
    sum_trpz = 0;
    for q = 1:N
        X = a + ((q-1)*dx); 
        Y = sqrt(((X^2)-(c^2)))/X;
        if (q == 0)||(q == N)
            coeff = 0.5;
        else
            coeff = 1;
        end
    sum_trpz = sum_trpz + (coeff*Y*dx);
    end
    Int_trpz(N) = sum_trpz;
    
    % Simpson's Rule
    sum_simp = 0;
    for r = 1:N
         X = a + ((r-1)*dx); 
         Y = sqrt(((X^2)-(c^2)))/X;
         if (mod(r,2) == 0)
            coeff = 4;
         else
            coeff = 2;
         end
         if (r == 0)||(r == N)
            coeff = 1;
         end
    sum_simp = sum_simp + (coeff*Y*dx*(1/3));
    end
    Int_simp(N) = sum_simp;
end

ExtInt = ((sqrt(b^2 - c^2)) - (abs(c)*asec(b/c))) - ((sqrt(a^2 - c^2)) - (abs(c)*asec(a/c)));
icf = ExtInt*ones(1,100);
plot(Int_mdp,'-ob');hold on;
plot(Int_trpz,'-r','linewidth',2);hold on;
plot(Int_simp,'-sm','linewidth',2);hold on;
plot(icf,'*k');hold on;
legend('Midpoint','Trapeziodal','Simpson','Real Value');
title('Implementation of Mid point, Trapeziodal and Simpson Rule','fontsize',14);
xlabel('Number of Interval->','fontsize',14);
ylabel ('Error Events ->','fontsize',14);
grid on;