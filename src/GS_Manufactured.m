clear; clc;
% The manufactured solution uses U=cos(Hx)*cos(Dy)
%Based on the problem statement, this gives F=(D^2+H^2)*cos(Hx)*cos(Dy)
%% Checkpoint Check
if exist( 'checkpointGS_Man.mat','file' ) % If a checkpoint file exists, load it
    load('checkpointGS_Man.mat')
end
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
W=zeros(Nx+2,Ny+2); %Preallocate a dummy matrix
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
%% SOR Loop
Count=0;    %Initialize the count
Max=1;  %Set Max greater than the limit
End=Nx+2;   %Precompute 'N'
while Max > 10^-8
    for k=2:Ny+1    %All y points not on the boundary
        y=Hy*(k-1)+ay;  %Compute the y-value for the given k
        for j=2:Nx+1    %All x points
            x=Hx*(j-1)+ax;  %Compute the x-value for the given j
            F=-(-D^2-H^2)*cos(D*x)*cos(H*y);  %Define F(x,y) for the particular j,k point
            U(j,k)=(Dx*(U(j,k-1)+U(j,k+1))+Dy*(U(j-1,k)+U(j+1,k))+(Dx*Dy*F))/(2*(Dx+Dy));   %Compute U at all interior points
        end
    end
    Count=Count+1;  %Increase the count
    Max=max(max(abs((W-U)./W)));  %Find the overall max
    W=U;    %n+1 becomes n
    if mod(Count,1000)==0   %Save checkpoint file every 1000 iterations
        save('checkpointGS_Man.mat'); %Save the file
    end     %Close if loop
end     %Close while loop
%% 3D Plot of the Matrix
X=-pi:Hx:pi;    %Discretize the X axis
Y=-pi:Hy:pi;    %Discretize the Y axis
V=transpose(U); %Transpose the matrix so that the x and y axes are correct
S=zeros(Nx+2,Ny+2); %Preallocate matrix S
for j=1:Ny+2
    for i=1:Nx+2
        S(i,j)=cos(H*X(i))*cos(D*Y(j)); %Determine the correct values for the manufactured solution
    end
end
T=transpose(S); %Transpose the S matrix for plotting
ERROR=zeros(Nx+2,Ny+2); %Preallocate the ERROR matrix
for j=2:Ny+1
    for i=2:Nx+1
        ERROR(i,j)=((U(i,j)-S(i,j)))^2; %Determine the squared difference at all interior points
    end
end
L2norm=sqrt(sum(sum(ERROR))/(Nx*Ny));   %Compute the L2 error
logL2=log10(L2norm);    %Determine log of the L2 norm for spacial accuracy
logdelta=log10(Hx); %Determine log of delta x for spacial accuracy
figure()    %First figure
h=surf(X,Y,V);  %Create surface plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines
figure()    %Second figure
zlevels=-0.8:0.4:0.8;   %Set the levels for the contour lines
contour(X,Y,V,zlevels,'ShowText','on');  %Create contour plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines
if exist( 'checkpointGS_Man.mat','file' )
    delete('checkpointGS_Man.mat');    %Delete checkpoint file once evertything is complete
end
