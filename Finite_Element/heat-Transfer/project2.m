%Boundary Element Method
clear all;close all;
N_x=4;
N=N_x;
h=2/N_x;
%coordinates of element endpoints
x_elements=[];
for j=1:N_x
    x_elements=[x_elements; -2 2-h*(j-1)];
end

for j=1:N_x
    x_elements=[x_elements; -2+h*(j-1) 0];
end

for j=1:N_x
    x_elements=[x_elements; 0 -h*(j-1)];
end

for j=1:N_x
    x_elements=[x_elements; h*(j-1) -2];
end

for j=1:2*N_x
    x_elements=[x_elements; 2 -2+h*(j-1)];
end

for j=1:2*N_x
    x_elements=[x_elements; 2-h*(j-1) 2];
end
x_elements=[x_elements; -2 2];

N_elements=8*N_x;

%coordinates of nodal points
epsilon=0.05;
epsilon0=epsilon*h;
x_nodes=[];
for j=1:(N_x+1)
    x_nodes=[x_nodes; -2 2-(j-1)*h];
end

%shifting
x_nodes(1,2)=x_nodes(1,2)-epsilon0;        

x_nodes(N_x+1,2)=x_nodes(N_x+1,2)+epsilon0;
for j=1:(N_x+1)
    x_nodes=[x_nodes; -2+(j-1)*h 0];       
end
x_nodes(N_x+2,1)=x_nodes(N_x+2,1)+epsilon0;
x_nodes(2*N_x+2,1)=x_nodes(2*N_x+2,1)-epsilon0; 
    x_nodes=[x_nodes; 0 -(j-1)*h];
end
x_nodes(2*N_x+3,2)=x_nodes(2*N_x+3,2)-epsilon0;
x_nodes(3*N_x+3,2)=x_nodes(3*N_x+3,2)+epsilon0;
for j=1:(N_x+1)
    x_nodes=[x_nodes; (j-1)*h -2];
end
x_nodes(3*N_x+4,1)=x_nodes(3*N_x+4,1)+epsilon0;
x_nodes(4*N_x+4,1)=x_nodes(4*N_x+4,1)-epsilon0;
for j=1:(2*N_x+1)
    x_nodes=[x_nodes; 2 -2+(j-1)*h];
end
x_nodes(4*N_x+5,2)=x_nodes(4*N_x+5,2)+epsilon0;
x_nodes(6*N_x+5,2)=x_nodes(6*N_x+5,2)-epsilon0;
for j=1:(2*N_x+1)
    x_nodes=[x_nodes; 2-(j-1)*h 2];
end
x_nodes(6*N_x+6,1)=x_nodes(6*N_x+6,1)-epsilon0;
x_nodes(8*N_x+6,1)=x_nodes(8*N_x+6,1)+epsilon0;
%plot(x_nodes(:,1),x_nodes(:,2),'o','MarkerSize',10);
N_nodes=8*(N_x+1)-2;

%the first and second node numbers and type of discontinuity
Element_nodes=[];
for j=1:N_x
    Element_nodes=[Element_nodes; j j+1 0];
end
Element_nodes(1,3)=1;
Element_nodes(N_x,3)=-1;
for j=1:N_x
    Element_nodes=[Element_nodes; j+N_x+1 j+N_x+2 0];
end
Element_nodes(N_x+1,3)=1;
Element_nodes(N_x+N_x,3)=-1;
for j=1:N_x
    Element_nodes=[Element_nodes; j+2*N_x+2 j+2*N_x+3 0];
end
Element_nodes(2*N_x+1,3)=1;
Element_nodes(2*N_x+N_x,3)=-1;
for j=1:N_x
    Element_nodes=[Element_nodes; j+3*N_x+3 j+3*N_x+4 0];
end
Element_nodes(3*N_x+1,3)=1;
Element_nodes(3*N_x+N_x,3)=-1;
for j=1:2*N_x
    Element_nodes=[Element_nodes; j+4*N_x+4 j+4*N_x+5 0];
end
Element_nodes(4*N_x+1,3)=1;
Element_nodes(5*N_x+N_x,3)=-1;
for j=1:2*N_x
    Element_nodes=[Element_nodes; j+6*N_x+5 j+6*N_x+6 0];
end
Element_nodes(6*N_x+1,3)=1;
Element_nodes(7*N_x+N_x,3)=-1;

%specify boundary conditions
Dirichlet_BC=[];
for i=1:N_x+1
    b_dirichlet(i)=2*(2-(i-1)*h)-(2-h*(i-1))^2;
end
for j=1:N_x+1
    Dirichlet_BC=[Dirichlet_BC; j b_dirichlet(j)];
end
N_Dirichlet=(N_x+1);
Newmann_BC=[];
N_Newmann=7*(N_x+1)-2;
for j=N_Dirichlet:N_nodes
    Newmann_BC=[Newmann_BC; j 0];
end

%------SOLVE THE BOUNDARY INTEGRAL EQUATION-------------------------

%form normal vectors
n_vector=[];
L=(1:N_elements)*0;
for j=1:N_elements
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    L(j)=norm(beta);
    n_vector=[n_vector; beta(2)/L(j) -beta(1)/L(j)];
end

u=zeros(N_nodes,1);
q=zeros(N_nodes,1);

for j=1:N_Dirichlet
    u(Dirichlet_BC(j,1))=Dirichlet_BC(j,2);
end
for j=1:N_Newmann
    q(Newmann_BC(j,1))=Newmann_BC(j,2);
end

%choose a quadrature
N_p=8;
w=[0.101228536290376 0.222381034453374 0.313706645877887 0.362683783378362 0.362683783378362 0.313706645877887 0.222381034453374 0.101228536290376];
p=[-0.960289856497536 -0.796666477413627 -0.525532409916329 -0.183434642495650 0.183434642495650 0.525532409916329 0.796666477413627 0.960289856497536];
w_new=w/2;
p_new=(p+1)/2;
tmp0=0*(1:N_p);
tmp1=0*(1:N_p);
Gij=[0 0];
Hij=[0 0];
%form G and H matrices
G=zeros(N_nodes);
H=zeros(N_nodes);
Q_plus=[1 -epsilon;
    0 1-epsilon]/(1-epsilon);
Q_minus=[1-epsilon 0;
    -epsilon 1]/(1-epsilon);
for j=1:N_elements
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    j_1=Element_nodes(j,1);
    j_2=Element_nodes(j,2);
    
    for i=1:N_nodes
        alpha=x_elements(j,1:2)-x_nodes(i,1:2);
        if(i == j_1)
            
            if(Element_nodes(j,3) == 1)
                Gij(1)=-(L(j)/2/pi)*(0.5*log(L(j))-3/4+0.5*((1-epsilon)^2*log(1-epsilon)+(2*epsilon-epsilon^2)*log(epsilon)+epsilon));
                Gij(2)=-(L(j)/2/pi)*(0.5*log(L(j))-1/4+0.5*((1-epsilon^2)*log(1-epsilon)+epsilon^2*log(epsilon)-epsilon));
                Hij(1)=0;
                Hij(2)=0;
            else
                Gij(1)=-(L(j)/2/pi)*(0.5*log(L(j))-3/4);
                Gij(2)=-(L(j)/2/pi)*(0.5*log(L(j))-1/4);
                Hij(1)=0;
                Hij(2)=0;
            end
        elseif(i == j_2)
            if(Element_nodes(j,3) == -1)
                Gij(1)=-(L(j)/2/pi)*(0.5*log(L(j))-1/4+0.5*((1-epsilon^2)*log(1-epsilon)+epsilon^2*log(epsilon)-epsilon));
                Gij(2)=-(L(j)/2/pi)*(0.5*log(L(j))-3/4+0.5*((1-epsilon)^2*log(1-epsilon)+(2*epsilon-epsilon^2)*log(epsilon)+epsilon));
                Hij(1)=0;
                Hij(2)=0;
            else
                Gij(1)=-(L(j)/2/pi)*(0.5*log(L(j))-1/4);
                Gij(2)=-(L(j)/2/pi)*(0.5*log(L(j))-3/4);
                Hij(1)=0;
                Hij(2)=0;
            end
        else
            for k=1:N_p
                tmp0(k)=log(norm(alpha+beta*p_new(k)))*(1-p_new(k));
                tmp1(k)=log(norm(alpha+beta*p_new(k)))*p_new(k);
            end
            Gij(1)=-(L(j)/2/pi)*(w_new*tmp0');
            Gij(2)=-(L(j)/2/pi)*(w_new*tmp1');
            
            for k=1:N_p
                tmp0(k)=(1/(norm(alpha+beta*p_new(k)))^2)*(1-p_new(k));
                tmp1(k)=(1/(norm(alpha+beta*p_new(k)))^2)*p_new(k);
            end
            Hij(1)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp0');
            Hij(2)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp1');
        end
        
        if(Element_nodes(j,3) == 1)
            Gij=Gij*Q_plus;
            Hij=Hij*Q_plus;
        end
        if(Element_nodes(j,3) == -1)
            Gij=Gij*Q_minus;
            Hij=Hij*Q_minus;
        end
        G(i,j_1)=G(i,j_1)+Gij(1);
        G(i,j_2)=G(i,j_2)+Gij(2);
        H(i,j_1)=H(i,j_1)+Hij(1);
        H(i,j_2)=H(i,j_2)+Hij(2);
    end
end

H=H+diag(0.5*ones(N_nodes,1));

for j=1:N_Dirichlet
    jd=Dirichlet_BC(j,1);
    q(jd,1)=u(jd,1);
    for i=1:N_nodes
        swap=H(i,jd);
        H(i,jd)=-G(i,jd);
        G(i,jd)=-swap;
    end
end
b=G*q;
u=H\b;
for j=1:N_Dirichlet
    jd=Dirichlet_BC(j,1);
    swap=u(jd,1);
    u(jd,1)=q(jd,1);
    q(jd,1)=swap;
end
length(u)
%---------POST PROCESSING-----------------------------------------------

%coordinates of nodal points

%inner nodal points
Nodes=[];
for j=1:N_x-1
    for i=1:2*N_x-1
        Nodes=[Nodes;i*h-2 2-j*h];
    end
end
for j=0:N_x-1
    for i=1:N_x-1
        Nodes=[Nodes;i*h -j*h];
    end
end
%plot(Nodes(:,1),Nodes(:,2),'o','MarkerSize',10);

Nodes=[Nodes;x_nodes]; %append all the boundary nodes

plot(Nodes(:,1),Nodes(:,2),'o','MarkerSize',10);

%choose the points for output of the solution
% x_output=Nodes(:,1);
% y_output=Nodes(:,2);     %for all nodal points
 x_output=0.5:0.5:1.5;
 y_output=0.5:0.5:1.5;
rowNum=size(Nodes,1);
u_output=zeros(1,3);      %or zeros(1,rowNum) if you wanna have solutions for all nodal points
for i=1:3                 %or i=1:rowNum      if you wanna have solutions for all nodal points
    xp=x_output(i);
    yp=y_output(i);
    g_vector=0*(1:N_nodes);
    h_vector=0*(1:N_nodes);
    for j=1:N_elements
        beta=x_elements(j+1,1:2)-x_elements(j,1:2);
        alpha=x_elements(j,1:2)-[xp yp];
        j_1=Element_nodes(j,1);
        j_2=Element_nodes(j,2);
        gj=[0 0];
        hj=[0 0];
        %g vector
        for k=1:N_p
            tmp0(k)=log(norm(alpha+beta*p_new(k)))*(1-p_new(k));
            tmp1(k)=log(norm(alpha+beta*p_new(k)))*p_new(k);
        end
        gj(1)=-(L(j)/2/pi)*(w_new*tmp0');
        gj(2)=-(L(j)/2/pi)*(w_new*tmp1');
        
        %H matrix
        for k=1:N_p
            tmp0(k)=1/(norm(alpha+beta*p_new(k)))^2*(1-p_new(k));
            tmp1(k)=1/(norm(alpha+beta*p_new(k)))^2*p_new(k);
        end
        hj(1)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp0');
        hj(2)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp1');
        if(Element_nodes(j,3) == 1)
            gj=gj*Q_plus;
            hj=hj*Q_plus;
        end
        if(Element_nodes(j,3) == -1)
            gj=gj*Q_minus;
            hj=hj*Q_minus;
        end
        g_vector(j_1)=g_vector(j_1)+gj(1);
        g_vector(j_2)=g_vector(j_2)+gj(2);
        h_vector(j_1)=h_vector(j_1)+hj(1);
        h_vector(j_2)=h_vector(j_2)+hj(2);
    end
    
    u_output(i)=g_vector*q-h_vector*u;
end
N_x
A_value=u_output(1)         % for the selected points
B_value=u_output(2)
C_value=u_output(3)

% x = Nodes(:,1);           % x,y,z for all points
% y = Nodes(:,2);
% z = u_output;


