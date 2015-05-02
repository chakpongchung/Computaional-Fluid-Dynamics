%------------INPUT GEOMETRY AND BOUNDRY CONDITIONS
clear;close all;
%coordinates of element endpoints
x_elements=[1 0; 1.25 0; 1.5 0; 1.75 0;
            2 0;2*cos(pi/8) 2*sin(pi/8); 2*cos(pi/4) 2*sin(pi/4); 2*cos(3*pi/8) 2*sin(3*pi/8);
            0 2; 0 1.75; 0 1.5; 0 1.25;
            0 1;cos(3*pi/8) sin(3*pi/8); cos(pi/4) sin(pi/4);cos(pi/8) sin(pi/8);
            1 0];
        plot(x_elements(:,1),x_elements(:,2),'o','MarkerSize',10);
N_nodes=16;

%specify boundary conditions
Dirichlet_BC=[5 20;
              6 20;
              7 20;
              8 20;
              13 10;
              14 10;
              15 10;
              16 10];
N_Dirichlet=8;
Newmann_BC=[1 0;
            2 0;
            3 0;
            4 0;
            9 0;
            10 0;
            11 0;
            12 0];
N_Newmann=8;

%------SOLVE THE BOUNDARY INTEGRAL EQUATION-------------------------

%coordinates of nodal points
x_nodes=[];
for j=1:N_nodes
    x_nodes=[x_nodes; 0.5*(x_elements(j,1:2)+x_elements(j+1,1:2))];
end

%form normal vectors
n_vector=[];
L=(1:N_nodes)*0;
for j=1:N_nodes
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    L(j)=norm(beta);
    n_vector=[n_vector; beta(2)/L(j) -beta(1)/L(j)];
end


u=zeros(N_nodes,1);
q=zeros(N_nodes,1);

for j=1:N_Dirichlet
    u(Dirichlet_BC(j,1),1)=Dirichlet_BC(j,2);
end
for j=1:N_Newmann
    q(Newmann_BC(j,1),1)=Newmann_BC(j,2);
end

%choose a quadrature
N_p=4;
w=[0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];
p=[-0.861136311594053 -0.339981043584856 0.339981043584856 0.861136311594053];
w_new=w/2;
p_new=(p+1)/2;
tmp2=0*(1:N_p);
%form G and H matrices
for j=1:N_nodes
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    for i=1:N_nodes
        if(i ~= j)
            alpha=x_elements(j,1:2)-x_nodes(i,1:2);
            
            %G matrix
            for k=1:N_p
                tmp2(k)=log(norm(alpha+beta*p_new(k)));
            end
            G(i,j)=-(L(j)/2/pi)*(w_new*tmp2');
            
            %H matrix
            for k=1:N_p
                tmp2(k)=1/(norm(alpha+beta*p_new(k)))^2;
            end
            H(i,j)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp2');
            
        else
            G(i,i)=-(L(j)/2/pi)*(-1+log(L(j)/2));
            H(i,i)=0.5;
        end
    end
end

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

%---------POST PROCESSING-----------------------------------------------

%choose the points for output of the solution
nx=4;
x_output=1.5*cos(pi/8*(-1/2+1:nx));
y_output=1.5*sin(pi/8*(-1/2+1:nx));
u_output=0*(1:nx);


for i=1:nx
    xp=x_output(i);
    yp=y_output(i);
    for j=1:N_nodes
        beta=x_elements(j+1,1:2)-x_elements(j,1:2);
        alpha=x_elements(j,1:2)-[xp yp];
        %g vector
        for k=1:N_p
            tmp2(k)=log(norm(alpha+beta*p_new(k)));
        end
        g_vector(j)=-(L(j)/2/pi)*(w_new*tmp2');
            
        %H matrix
        for k=1:N_p
            tmp2(k)=1/(norm(alpha+beta*p_new(k)))^2;
        end
        h_vector(j)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp2');
    end
    u_output(i)=g_vector*q-h_vector*u;
end
u_output


%for contour plots
vortices=[];
for i=1:N_nodes
    vortices=[vortices; x_elements(i,1:2)];
end
vortices=[vortices;1.5*cos(pi/16) 1.5*sin(pi/16); 1.5*cos(3*pi/16) 1.5*sin(3*pi/16); 
    1.5*cos(5*pi/16) 1.5*sin(5*pi/16);1.5*cos(7*pi/16) 1.5*sin(7*pi/16)];
B=[1 2 16;
    2 16 17;
    2 3 17;
    3 4 17;
    4 6 17;
    4 5 6;
    15 16 18;
    16 17 18;
    6 17 18;
    6 7 18;
    14 15 19;
    15 18 19;
    7 18 19;
    7 8 19;
    14 19 20;
    8 19 20;
    13 12 14;
    12 14 20;
    12 11 20;
    11 10 20;
    10 8 20;
    10 9 8];
   
u_output=u;

for i=N_nodes+1:N_nodes+4
    xp=vortices(i,1);
    yp=vortices(i,2);
    for j=1:N_nodes
        beta=x_elements(j+1,1:2)-x_elements(j,1:2);
        alpha=x_elements(j,1:2)-[xp yp];
        %g vector
        for k=1:N_p
            tmp2(k)=log(norm(alpha+beta*p_new(k)));
        end
        g_vector(j)=-(L(j)/2/pi)*(w_new*tmp2');
            
        %H matrix
        for k=1:N_p
            tmp2(k)=1/(norm(alpha+beta*p_new(k)))^2;
        end
        h_vector(j)=-(L(j)/2/pi)*(alpha*n_vector(j,1:2)')*(w_new*tmp2');
    end
    u_output=[u_output; g_vector*q-h_vector*u];
end

%figure(2)
patch('Vertices',vortices,'Faces',B,'FaceVertexCData',u_output,'FaceColor','interp');