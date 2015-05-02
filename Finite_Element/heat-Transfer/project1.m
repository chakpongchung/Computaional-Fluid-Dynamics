%Finite Element Method
clear all; close all;
N_x=4;
h=2/N_x;
%N_e: number of elements
%N_p: number of nodes, including boundary nodes
%N_b: number of boundary nodes of Dirichlet type
N_e=6*N_x*N_x;
N_p=(N_x*2+1)^2-N_x^2;
N_b=N_x+1;
b_dirichlet=zeros(N_x+1,1);
for i=1:N_x+1
    b_dirichlet(i)=2*(2-(i-1)*h)-(2-h*(i-1))^2;
end

%coordinates of nodal points
Nodes=[];

for j=0:N_x
    for i=0:2*N_x
        Nodes=[Nodes;i*h-2 2-j*h];
    end
end
for j=1:N_x
    for i=0:N_x
        Nodes=[Nodes;i*h -j*h];
    end
end

 plot(Nodes(:,1),Nodes(:,2),'o','MarkerSize',10);
 
 %B is the connectivity matrix
B=[];
for j=1:N_x
    for i=1:2*N_x
        B=[B;(j-1)*(2*N_x+1)+i (j-1)*(2*N_x+1)+i+1 j*(2*N_x+1)+i+1];
        B=[B;j*(2*N_x+1)+i+1 j*(2*N_x+1)+i (j-1)*(2*N_x+1)+i];
    end
end
for j=1:N_x
    for i=1:N_x
        B=[B;(j-1)*(N_x+1)+i+N_x*(2*N_x+2) (j-1)*(N_x+1)+i+1+N_x*(2*N_x+2) (j-1)*(N_x+1)+i+2+N_x*(2*N_x+3)];
        B=[B;(j-1)*(N_x+1)+i+2+N_x*(2*N_x+3) (j-1)*(N_x+1)+i+1+N_x*(2*N_x+3) (j-1)*(N_x+1)+i+N_x*(2*N_x+2)];
    end
end

Dirichlet_nodes=[];
for j=0:N_x
    Dirichlet_nodes=[Dirichlet_nodes; j*(2*N_x+1)+1];
end
A_element=[1/2 -1/2 0;
          -1/2 1 -1/2; 
           0 -1/2 1/2];

%----assemble the global matrix and solution is stored in c_global----------------
A_global=zeros(N_p);
b_global=zeros(N_p,1);

for i_e=1:N_e
    for m=1:3
        for n=1:3
            i=B(i_e,m);
            j=B(i_e,n);
            A_global(i,j)=A_global(i,j)+A_element(m,n);
        end
    end
end

b_global(Dirichlet_nodes)=b_dirichlet;
A_global(Dirichlet_nodes,:)=0;
for i=1:(N_x+1)
    A_global(Dirichlet_nodes(i),Dirichlet_nodes(i))=1;
end
c_global=A_global\b_global;

%Solutions at indicated points:
A_node=3/4*N_x*(2*N_x+1)+5/4*N_x+1;
B_node=2/4*N_x*(2*N_x+1)+3/2*N_x+1;
C_node=1/4*N_x*(2*N_x+1)+7/4*N_x+1;
N_x
A_value=c_global(A_node)
B_value=c_global(B_node)
C_value=c_global(C_node)


