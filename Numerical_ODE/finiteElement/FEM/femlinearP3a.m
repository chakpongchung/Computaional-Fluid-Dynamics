clear all

%Godfred
nvals = [9 19 29];
%end Godfred

ntest = length(nvals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);
for jtest=1:ntest
    N=nvals(jtest);
    hvals(jtest)=1/nvals(jtest);
    h=hvals(jtest);
    x=(1:N)'*h;
    
    %f=x.*(1-x)
    %   GODFRED f
    f=-(cos(pi*x)).^2 - 2*pi^2*cos(2*pi*x);
    % end  Godfred f
    
    beta=0;
    e=exp(1);
    
    %exactu=3/2*exp(x-1)+2*exp(-x)-3/2*exp(-x-1)-x.^2+x-2
    %  Godfred exactu
    exactu = exp(x)/(1+e) + e*exp(-x)/(1+e) - (cos(pi*x)).^2;
    %end Godfred
    
    A=zeros(N,N);
    for i=1:N-1
        A(i,i)=2/h+2*h/3;
        A(i,i+1)=-1/h+h/6;
        A(i+1,i)=-1/h+h/6;
    end
    A(N,N)=1/h+h/3+1;
    
    %  Godfred
    %A(1,1) = 0;
    A(N,N) = 0;
    %end Godfred
    
    b=h*f;
    b(N)=b(N)/2+beta;
    
    %  Godfred
    %b(1) = 0;
    b(N) = 0;
    % end
    
    
    c=A\b;
    [0;c];
    for i=1:N
        L2(i)=(exactu(i)-c(i)).^2;
    end
    L2norm(1)=0;
    for i=1:(N-1)
        L2norm(i+1)=L2norm(1)+L2(i+1);
    end
    E(jtest)=sqrt(L2norm(N));
    figure(jtest)
    plot([0;x],[0;c],'--o',[0;x],[0;exactu],'r-')
    % hold on
    %legend('FEM','exact',1);
    
    %Godfred comment out
    %axis([0 1 0 0.1])
    
    xlabel('x')
    ylabel('u')
    title(['N=',num2str(N),', L2 norm ',num2str(sqrt(h)*norm(c-exactu,2))])
    % hold on
end
%figure(4)
%error_table
[hvals, E]
%error_log
[hvals, log(E)]