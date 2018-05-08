clear; clc;
%% Checkpoint Check
if exist( 'checkpointSOR_zero.mat','file' ) % If a checkpoint file exists, load it
    load('checkpointSOR_zero.mat')
end
%% Define basic parameters
ax=-pi; %Define lower x bound
ay=-pi; %Define lower y bound
bx=pi;  %Define upper x bound
by=pi;  %Define upper y bound
Lx=bx-ax;  %Define length of x
Ly=by-ay;  %Define length of y
%% Initial Matrix
Nx=160;  %Number of nodes added to the x-axis
Ny=160;  %Number of nodes added to the y-axis
Hx=Lx/(1+Nx);    %Length of x-axis segment
Dx=Hx*Hx;   %Determine delta x squared
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
%% SOR Loop
Count=0;    %Initialize the count
Max=1;  %Set Max greater than the limit
w=1.9;  %Set the coefficient for the over-relaxation
End=Nx+2;   %Precompute 'N'
while Max > 10^-6
    for k=2:Ny+1    %All y points not on the boundary
        y=Hy*(k-1)+ay;  %Compute the y-value for the given k
        %BC for x=ax
        U(1,k)=(Dx*(U(1,k-1)+U(1,k+1))+Dy*2*U(2,k))/(2*(Dx+Dy));  %Compute boundary value at x=ax
        %BC for x=bx
        F=0;  %Define F(x,y) for the particular j,k point
        U(End,k)=((Dx*(U(End,k-1)+U(End,k+1))+Dy*2*U(End-1,k))+(Dx*Dy*F))/(2*(Dx+Dy)); %Compute boundary value at x=bx
        for j=2:Nx+1    %All x points
            x=Hx*(j-1)+ax;  %Compute the x-value for the given j
            F=0;  %Define F(x,y) for the particular j,k point
            U(j,k)=-(w-1)*W(j,k)+w*(Dx*(U(j,k-1)+U(j,k+1))+Dy*(U(j-1,k)+U(j+1,k))+(Dx*Dy*F))/(2*(Dx+Dy));
        end
    end
    AverageValue=sum(sum(U))/((Nx+2)*(Ny+2));
    Count=Count+1;  %Increase the count
    Max=max(max(abs((W-U)./W)));  %Find the overall max
    W=U;    %n+1 becomes n
    if mod(Count,1000)==0   %Save checkpoint file every 1000 iterations
        save('checkpointSOR_zero.mat'); %Save the file
    end     %Close if loop
end     %Close while loop for Max
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
zlevels=[-200,-150,-100,-50,-40,-30,-20,-10,0,10,20,30,40];
contour(X,Y,V,zlevels,'ShowText','on');  %Create contour plot
ylabel('y') %Label the y-axis
xlabel('x') %Label the x-axis
set(h,'linestyle','none');  %Remove the gridlines
if exist( 'checkpointSOR_zero.mat','file' )
    delete('checkpointSOR_zero.mat');    %Delete checkpoint file once evertything is complete
end