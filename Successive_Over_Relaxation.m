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
W=zeros(Nx+2,Ny+2);
for j=1:Nx+2
    x=Hx*(j-1)+ax;  %Compute the x-value for the given i
    U(j,1)=((x-ax)^2)*sin(pi/2*(x-ax)/(Lx));  %Set BC for y=ay
    U(j,Ny+2)=cos(pi*(x-ax))*cosh(bx-x);    %Set BC for y=by
end
%% SOR Loop
Count=0;
Max=1;
w=1.9;  %Set the coefficient for the over-relaxation
End=Nx+2;
while Max > 10^-6
    for k=2:Ny+1    %All y points not on the boundary
        y=Hy*(k-1)+ay;  %Compute the y-value for the given k
        %BC for x=ax
        U(1,k)=(Dx*(U(1,k-1)+U(1,k+1))+Dy*2*U(2,k))/(2*(Dx+Dy));  %Compute boundary value at x=ax
        %BC for x=bx
        F=sin(pi*(bx-ax)/Lx)*cos((0.5*pi)*(2*(y-ay)/Ly)+1);  %Define F(x,y) for the particular j,k point
        U(End,k)=((Dx*(U(End,k-1)+U(End,k+1))+Dy*2*U(End-1,k))+(Dx*Dy*F))/(2*(Dx+Dy)); %Compute boundary value at x=bx
        for j=2:Nx+1    %All x points
            x=Hx*(j-1)+ax;  %Compute the x-value for the given j
            F=sin(pi*(x-ax)/Lx)*cos((0.5*pi)*(2*(y-ay)/Ly)+1);  %Define F(x,y) for the particular j,k point
            U(j,k)=-(w-1)*W(j,k)+w*(Dx*(U(j,k-1)+U(j,k+1))+Dy*(U(j-1,k)+U(j+1,k))+(Dx*Dy*F))/(2*(Dx+Dy));
        end
    end
    Count=Count+1;  %ERROR
    Diff=W-U;
    MaxA=max(abs(Diff));
    Max=max(MaxA);
    W=U;
end
%% 3D Plot of the Matrix
X=-pi:Hx:pi;
Y=-pi:Hy:pi;
V=transpose(U);
h=surf(X,Y,V);
ylabel('y')
xlabel('x')
set(h,'linestyle','none');