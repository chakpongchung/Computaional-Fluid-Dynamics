options = odeset('RelTol',1e-5,'AbsTol',1e-5);
[T,Y] = ode45(@curtiss,[0 10],1,options);
y=2500/2501*cos(T)+50/2501*sin(T)+1/2501 * exp(-50*T);
plot(T,Y(:,1),'-',T(1:5:end),y(1:5:end),'o');

