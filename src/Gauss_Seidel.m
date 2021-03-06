clear; clc;
%% Checkpoint Check
if exist( 'checkpointGS.mat','file' ) % If a checkpoint file exists, load it
    load('checkpointGS.mat')
end
%% Define basic parameters
ax=-pi; %Define lower x bound
ay=-pi; %Define lower y bound
bx=pi;  %Define upper x bound
by=pi;  %Define upper y bound
Lx=bx-ax;  %Define length of X
Ly=by-ay;  %Define length of Y
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
for j=1:Nx+2
    x=Hx*(j-1)+ax;  %Compute the x-value for the given i
    U(j,1)=((x-ax)^2)*sin(pi/2*(x-ax)/(Lx));  %Set BC for y=ay
    U(j,Ny+2)=(cos(pi*(x-ax))-1)*cosh(bx-x);    %Set BC for y=by
end
%% Gauss-Seidel Loop
Count=0;    %Initialize the count
Max=1;  %Set Max greater than the limit
End=Nx+2;   %Precompute 'N'
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
            U(j,k)=(Dx*(U(j,k-1)+U(j,k+1))+Dy*(U(j-1,k)+U(j+1,k))+(Dx*Dy*F))/(2*(Dx+Dy));
        end
    end
    AverageValue=sum(sum(U))/((Nx+2)*(Ny+2));   %Calculate the average of U value
    Count=Count+1;  %Increase the count
    Max=max(max(abs((W-U)./W)));  %Find the overall max
    W=U;    %n+1 becomes n
    if mod(Count,1000)==0   %Save checkpoint file every 1000 iterations
        save('checkpointGS.mat'); %Save the file
    end     %Close if loop
end     %Close while loop
%% 3D Plot of the Matrix
X=-pi:Hx:pi;    %Discretize the X axis
Y=-pi:Hy:pi;    %Discretize the Y axis
V=transpose(U); %Transpose the matrix so that the x and y axes are correct
figure()    %First figure
h=surf(X,Y,V);  %Create surface plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines
figure()    %Second figure
zlevels=[-100,-80,-70,-50,-20,-10,-5,-2,0,10,20,30,40];    %Set the levels for the contour plot
contour(X,Y,V,zlevels,'ShowText','on');  %Create contour plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines
if exist( 'checkpointGS.mat','file' )
    delete('checkpointGS.mat');    %Delete checkpoint file once everything is complete
end
