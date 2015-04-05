function value=trap(f,a,b,n)
%TRAP: Approximates the integral of f from a to b
%using the Trapezoidal Rule with n subintervals
8%of equal width.
value = (f(a)+f(b))/2;
dx = (b-a)/n;
for k=1:(n-1)
c = a+k*dx;
value = value + f(c);
end
value = dx*value;

return