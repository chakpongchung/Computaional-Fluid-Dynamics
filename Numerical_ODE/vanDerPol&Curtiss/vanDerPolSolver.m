options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5]);
[T,Y] = ode45(@vanDerPol,[0 25],[0.5,0.5],options);
plot(T,Y(:,1),'-',T,Y(:,2),'-.')