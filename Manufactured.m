clear; clc;
% The manufactured solution uses U=cos(Hx)*cos(Dy)
%Based on the problem statement, this gives F=(D^2+H^2)*cos(Hx)*cos(Dy)
%% Define basic parameters
ax=-pi; %Define lower x bound
ay=-pi; %Define lower y bound
bx=pi;  %Define upper x bound
by=pi;  %Define upper y bound
Lx=bx-ax;  %Define length of X
Ly=by-ay;  %Define length of Y
H=2;    %Define the H value for the manufactured solution
D=2;    %Define the H value for the manufactured solution
%% Nodes
Nx=160;    %Number of nodes added to the x-axis
Hx=Lx/(1+Nx);    %Length of x-axis segment
Dx=Hx*Hx;   %Determine delta x squared
Ny=160;    %Number of nodes added to the y-axis
Hy=Ly/(1+Ny);    %Length of y-axis segment
Dy=Hy*Hy;   %Determine delta y squared
%% Create U matrix
U=zeros(Nx+2,Ny+2); %Preallocate the U matrix
for k=1:Ny+2
    y=Hy*(k-1)+ay;  %Compute the y-value for the given k
    for j=1:Nx+2
        x=Hx*(j-1)+ax;  %Compute the x-value for the given j
        U(j,1)=cos(H*x)*cos(D*y);  %Set BC for y=ay
        U(j,Ny+2)=cos(H*x)*cos(D*y);    %Set BC for y=by
        U(1,k)=cos(H*x)*cos(D*y);   %Set BC for x=ax
        U(Nx+2,k)=cos(H*x)*cos(D*y);    %Set BC for x=bx
    end
end
%% 3D Plot of the Matrix
X=-pi:Hx:pi;    %Discretize the X axis
Y=-pi:Hy:pi;    %Discretize the Y axis
S=zeros(Nx+2,Ny+2); %Preallocate matrix S
for j=1:Ny+2
    for i=1:Nx+2
        S(i,j)=cos(H*X(i))*cos(D*Y(j)); %Determine the correct values for the manufactured solution
    end
end
T=transpose(S); %Transpose the S matrix for plotting
figure()    %First figure
h=surf(X,Y,T);  %Create surface plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines
figure()    %Second figure
zlevels=-0.8:0.4:0.8;   %Set the levels for the contour lines
contour(X,Y,T,zlevels,'ShowText','on');  %Create contour plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines