clear;clc

for num = 1:2
    figure(num)
    J = [20;40;80;200];
    for k = 1:4
        theta = 0.6;
        N = 1000;
        T = 1;
        dt = T/N;
        t = 0:dt:T;
        X = 1;
        dx = X/J(k);
        x = 0:dx:1;
        mu = dt/(dx^2);
        U = zeros(J(k)+1,N+1);
        switch num
            case 1
                u0 = @(x) (x>=0.4).*(x<=0.6);
            case 2
                u0 = @(x) exp(-20*(x-1/2).^2);
        end
        U(2:1:J(k),1) = u0(dx:dx:(1-dx));
        
        for n=2:(N+1)
            temp = zeros(J(k)-1,1);
            for i = 1:(J(k)-1)
                temp(i) = mu*(1-theta)*U(i,n-1) + (1-2*mu*(1-theta))*U(i+1,n-1) + mu*(1-theta)*U(i+2,n-1);
            end
            U(2:1:J(k),n) = tridiag_solve(repmat(-mu*theta,J(k)-2,1),repmat(1+2*mu*theta,J(k)-1,1),repmat(-mu*theta,J(k)-2,1),temp);
        end
        

        subplot(2,2,k)
        hold on;
        title(sprintf('J=%d',J(k)))
        for i=1:20:(N+1)
            if i==1
                plot(x, U(:,i),'-b');
            else
                plot(x, U(:,i),'--r');
            end            
        end
        s1 = '(x>=0.4)(x<=0.6)';
        s2 = 'exp^{20(x-1/2)^2}';
        
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

        switch num
            case 1
                text(0.5, 1,sprintf('u = %s',s1),'HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontWeight','bold')
            case 2
                text(0.5, 1,sprintf('u = %s',s2),'HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontWeight','bold')
        end        
    end
end