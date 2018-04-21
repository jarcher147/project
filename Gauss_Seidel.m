clear; clc;
%% Define basic parameters
ax=-pi; %Define lower x bound
ay=-pi; %Define lower y bound
bx=pi;  %Define upper x bound
by=pi;  %Define upper y bound
Lx=bx-ax;  %Define length of X
Ly=by-ay;  %Define length of Y
%% Nodes
Nx=20;    %Number of nodes added to the x-axis
Hx=Lx/(1+Nx);    %Length of x-axis segment
Ny=20;    %Number of nodes added to the y-axis
Hy=Ly/(1+Ny);    %Length of y-axis segment
%% Create U matrix
U=zeros(Nx+2,Ny+2); %Preallocate the U matrix
for i=1:Nx+2
    x=Hx*(i-1)+ax;  %Compute the x-value for the given i
    U(1,i)=(x-ax)^2*sin(pi*(x-ax)/(2*Hx));  %Set BC for y=ay
    U(Nx+2,i)=cos(pi*(x-ax))*cosh(bx-x);    %Set BC for y=by
end
%% Gauss-Seidel Loop
while Count < 10000
    for j=2:Ny+1    %All y points not on the boundary
        y=Hy*(j-1)+ay;  %Compute the y-value for the given j
        
        
        for i=1:Nx+2    %All x points
            x=Hx*(i-1)+ax;  %Compute the x-value for the given i
            F=sin(pi*(x-ax)/Hx)*cos((0.5*pi)*(2*(y-ay)/Hy)+1);  %Define F(x,y) for the particular i,j point
            U(i,j)=(0.25)*(Hx*(U(i,j-1)+U(i,j+1))+Hy*(U(i-1,j)+U(i+1,j))+(Hx*Hy*F));
        end
    end
    Count=Count+1;  %ERROR
end