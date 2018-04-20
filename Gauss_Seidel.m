clear; clc;
%% Define basic parameters
Lx=pi;  %Define length of X
Ly=pi;  %Define length of Y
U0=0;   %Define boundary conditions
UL=0;
V0=0;
VL=0;
M=1;    %Integer in f(x)
%% Nodes
Nx=39;    %Number of nodes added to the x-axis
Hx=Lx/(1+Nx);    %Length of x-axis segment
Ny=39;    %Number of nodes added to the x-axis
Hy=Ly/(1+Ny);    %Length of y-axis segment
%% Create U matrix
U=zeros(Ny+2,Nx+2); %Create the initial U matrix
U(1,:)=U0;  %Apply boundary condition to bottom
U(Ny+2,:)=UL;  %Apply boundary condition to top
U(:,1)=V0;  %Apply boundary condition to left
U(:,Nx+2)=VL;  %Apply boundary condition to right
%% Gauss-Seidel Loop
Count=1;
while Count < 10000
    for k=Ny+1:-1:2    %All y points not on the boundary
        for j=2:Nx+1    %All x points not on the boundary
            F=-2*M*sin(M*(j-1)*Hx)*cosh(M*(Ny+2-k)*Hy);  %Define f(x,y) for the particular j,k point
            U(k,j)=(1/4)*(U(k,j-1)+U(k,j+1)+U(k-1,j)+U(k+1,j))-((Hx^2)/4*F);
        end
    end
    Count=Count+1;
end
%% Solve for the true value of U(x,y)
Utrue=zeros(Ny+2,Nx+2); %Create the initial U matrix
Utrue(1,:)=U0;  %Apply boundary condition to bottom
Utrue(Ny+2,:)=UL;  %Apply boundary condition to top
Utrue(:,1)=V0;  %Apply boundary condition to left
Utrue(:,Nx+2)=VL;  %Apply boundary condition to right
for y=Ny+1:-1:2
    for x=2:Nx+1
        Utrue(y,x)=(Ly-((Ny+2-y)*Hy))*sin(M*(x-1)*Hx)*sinh(M*(Ny+2-y)*Hy);
    end
end
Uerrorabs=abs(U(2:Ny+1,2:Nx+1)-Utrue(2:Ny+1,2:Nx+1));
Uerrorrel=abs(U(2:Ny+1,2:Nx+1)-Utrue(2:Ny+1,2:Nx+1))./abs(Utrue(2:Ny+1,2:Nx+1));