%------------INPUT GEOMETRY AND BOUNDRY CONDITIONS
clear all;close all;
%coordinates of element endpoints
x_elements=[1 0; 1.25 0; 1.5 0; 1.75 0;2 0;
            2*cos(pi/8) 2*sin(pi/8); 2*cos(pi/4) 2*sin(pi/4); 2*cos(3*pi/8) 2*sin(3*pi/8);
            0 2; 0 1.75; 0 1.5; 0 1.25;0 1;
            cos(3*pi/8) sin(3*pi/8); cos(pi/4) sin(pi/4);cos(pi/8) sin(pi/8);
            1 0];
N_elements=16;
%coordinates of nodal points
epsilon=0.05;
epsilon0=epsilon/4  %the u_output is sensitive to this epsilon.
x_nodes=[1+epsilon0 0; 1.25 0; 1.5 0; 1.75 0; 2-epsilon0 0;
         2 epsilon0; 2*cos(pi/8) 2*sin(pi/8); 2*cos(pi/4) 2*sin(pi/4); 2*cos(3*pi/8) 2*sin(3*pi/8); epsilon0 2;
         0 2-epsilon0; 0 1.75; 0 1.5; 0 1.25; 0 1+epsilon0;
         epsilon0 1; cos(3*pi/8) sin(3*pi/8); cos(pi/4) sin(pi/4);cos(pi/8) sin(pi/8); 1 epsilon0];
plot(x_nodes(:,1),x_nodes(:,2),'o','MarkerSize',10);
N_nodes=20;

%the first and second node numbers and type of discontinuity
Element_nodes=[1 2 1;
               2 3 0;
               3 4 0;
               4 5 -1;
               
               6 7 1;
               7 8 0;
               8 9 0;
               9 10 -1;
               
               11 12 1;
               12 13 0;
               13 14 0;
               14 15 -1;
               
               16 17 1;
               17 18 0;
               18 19 0;
               19 20 -1];


%plot(x_nodes(:,1),x_nodes(:,2),'o','MarkerSize',10);

%specify boundary conditions
Dirichlet_BC=[6 20;
              7 20;
              8 20;
              9 20;
              10 20
              16 10;
              17 10;
              18 10;
              19 10;
              20 10];
N_Dirichlet=10;
Newmann_BC=[1 0;
            2 0;
            3 0;
            4 0
            5 0;
            
            11 0;
            12 0;
            13 0;
            14 0;
            15 0];
N_Newmann=10;

%------SOLVE THE BOUNDARY INTEGRAL EQUATION-------------------------

%form normal vectors
n_vector=[];
L=(1:N_elements)*0;
for j=1:N_elements
    beta=x_elements(j+1,1:2)-x_elements(j,1:2);
    L(j)=norm(beta);
    n_vector=[n_vector; beta(2)/L(j) -beta(1)/L(j)];
end
n_vector;

u=zeros(N_nodes,1);
q=zeros(N_nodes,1);

for j=1:N_Dirichlet
    u(Dirichlet_BC(j,1),1)=Dirichlet_BC(j,2);
end
for j=1:N_Newmann
    q(Newmann_BC(j,1),1)=Newmann_BC(j,2);
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
%---------POST PROCESSING-----------------------------------------------

%choose the points for output of the solution
N_x=4;
x_output=1.5*cos(pi/8*(-1/2+1:N_x));
y_output=1.5*sin(pi/8*(-1/2+1:N_x));%1/16,3/16,5/16,7/16
u_output=0*(1:N_x);
uuu=u;%% uuu is to agregate the output value for plotting
for i=1:N_x
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
    uuu=[uuu; g_vector*q-h_vector*u];%append the value u at each node
end
display('u_output');u_output
size(u_output)

Nodes=[];
for j=1:N_x
    for i=1:N_x
        Nodes=[Nodes;x_output(i) y_output(j)];
    end
end

%for contour plots
vortices=[];
for i=1:N_nodes
    vortices=[vortices; x_nodes(i,1:2)];
end
vortices=[vortices;1.5*cos(pi/16) 1.5*sin(pi/16); 1.5*cos(3*pi/16) 1.5*sin(3*pi/16); 
    1.5*cos(5*pi/16) 1.5*sin(5*pi/16);1.5*cos(7*pi/16) 1.5*sin(7*pi/16)];

%B is the connectivity matrix/ adjacency matrix to be determined
 B=[];
for j=1:N_x
    for i=1:N_x
        B=[B;(j)*(N_x+1)+i (j-1)*(N_x+1)+i (j)*(N_x+1)+i+1];
        B=[B;(j-1)*(N_x+1)+i+1 j*(N_x+1)+i+1 (j-1)*(N_x+1)+i];
    end
end  
patch('Vertices',vortices,'Faces',B,'FaceVertexCData',uuu,'FaceColor','interp');