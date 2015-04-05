close all;clear;clc;
T=sqrt(1/3);
h=T/2;
%N_e: number of elements
%N_p: number of nodes, including boundary nodes
%N_b: number of boundary nodes of Dirichlet type
N_e=16;
N_p=15;
N_b=12;
Nodes=[-T 0;
    -T/2 0;
    0 0;
    T/2 0;
    T 0;
    
    -3/4*T 0.25;
    -1/4*T 0.25;
    1/4*T 0.25;
    3/4*T 0.25;
    
    -T/2 0.5;
    0 0.5;
    T/2 0.5;
    
    -1/4*T 0.75;
    1/4*T 0.75;
    
    0 1];

%B is the connectivity matrix
B=[ 1 2 6;
    7 6 2;
    2 3 7;
    8 7 3;
    3 4 8;
    9 8 4;
    4 5 9;
    6 7 10;
    11 10 7;
    7 8 11;
    12 11 8;
    8 9 12;
    10 11 13;
    14 13 11;
    11 12 14;
    13 14 15];

% B=[  1 2 6;
%      2 6 7;
%      2 3 7;
%      3 7 8;
%      3 4 8;
%      4 8 9;
%      4 5 9;
%      7 10 6;
%      10 7 11;
%      8 11 7;
%      11 8 12;
%      9 12 8;
%      11 13 10;
%      13 11 14;
%      12 14 11;
%      14 15 13];

Dirichlet_nodes=[1
    2
    3
    4
    5
    6
    9
    10
    12
    13
    14
    15];

A_element=[T -1/2*T -1/2*T;
    -1/2*T T -1/2*T;
    -1/2*T  -1/2*T T];
% b_element=[h^2/6 h^2/6 h^2/6];
b_element=[h^2/4*T h^2/4*T h^2/4*T];


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
    for m=1:3
        i=B(i_e,m);
        b_global(i)=b_global(i)+b_element(m);
    end
end

map=(1:(N_p-N_b))';
i_map=1;
i_dirichlet=1;
for i=1:N_p
    if(i==Dirichlet_nodes(i_dirichlet))
        i_dirichlet=i_dirichlet+1;
    else
        map(i_map)=i;
        i_map=i_map+1;
    end
end
A_sub=zeros(N_p-N_b);
b_sub=zeros(N_p-N_b,1);
for i=1:N_p-N_b
    for j=1:N_p-N_b
        i0=map(i);
        j0=map(j);
        A_sub(i,j)=A_global(i0,j0);
    end
end
for i=1:N_p-N_b
    i0=map(i);
    b_sub(i,1)=b_global(i0,1);
end
c_sub=A_sub\b_sub;
c_sub
c_global=zeros(N_p,1);
for i=1:N_p-N_b
    i0=map(i);
    c_global(i0,1)=c_sub(i,1);
end
%----------------------------------------------------------------------------

%postprocessing
%plot the value at nodes 1 to 5
% exact=ones(5,1)*0.5;
% for n=1:20
%     alpha=(2*n-1)*pi/2;
%     exact=exact+2*(-1)^n*cosh(alpha*Nodes(1:5,1))/(alpha^3)/cosh(alpha);
% end
% plot(Nodes(1:5,1),c_global(1:5,1),'-o',Nodes(1:5,1),exact,'-s');
% xlabel('x');
% ylabel('u');
% legend('FEM','exact');
%figure(2)
%plot the distribution on each element
patch('Vertices',Nodes,'Faces',B,'FaceVertexCData',c_global,'FaceColor','interp');