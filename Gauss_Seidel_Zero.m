clear; clc;
%% Define basic parameters
ax=-pi; %Define lower x bound
ay=-pi; %Define lower y bound
bx=pi;  %Define upper x bound
by=pi;  %Define upper y bound
Lx=bx-ax;  %Define length of X
Ly=by-ay;  %Define length of Y
%% Nodes
Nx=80;    %Number of nodes added to the x-axis
Hx=Lx/(1+Nx);    %Length of x-axis segment
Dx=Hx*Hx;   %Determine delta x squared
Ny=80;    %Number of nodes added to the y-axis
Hy=Ly/(1+Ny);    %Length of y-axis segment
Dy=Hy*Hy;   %Determine delta y squared
%% Create U matrix
U=zeros(Nx+2,Ny+2); %Preallocate the U matrix
for j=1:Nx+2
    x=Hx*(j-1)+ax;  %Compute the x-value for the given i
    U(j,1)=((x-ax)^2)*sin(pi/2*(x-ax)/(Lx));  %Set BC for y=ay
    U(j,Ny+2)=cos(pi*(x-ax))*cosh(bx-x);    %Set BC for y=by
end
%% Gauss-Seidel Loop
Count=0;
End=Nx+2;
while Count < 10000
    for k=2:Ny+1    %All y points not on the boundary
        y=Hy*(k-1)+ay;  %Compute the y-value for the given k
        %BC for x=ax
        U(1,k)=(Dx*(U(1,k-1)+U(1,k+1))+Dy*2*U(2,k))/(2*(Dx+Dy));  %Compute boundary value at x=ax
        %BC for x=bx
        
        U(End,k)=(Dx*(U(End,k-1)+U(End,k+1))+Dy*2*U(End-1,k))/(2*(Dx+Dy)); %Compute boundary value at x=bx
        for j=2:Nx+1    %All x points
            x=Hx*(j-1)+ax;  %Compute the x-value for the given j
            
            U(j,k)=(Dx*(U(j,k-1)+U(j,k+1))+Dy*(U(j-1,k)+U(j+1,k)))/(2*(Dx+Dy));
        end
    end
    Count=Count+1;  %ERROR
end
%% 3D Plot of the Matrix
X=-pi:Hx:pi;
Y=-pi:Hy:pi;
V=transpose(U);
h=surf(X,Y,V);
ylabel('y')
xlabel('x')
set(h,'linestyle','none');
colormap('gray');