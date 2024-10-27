% To compute stress formulation with Dirichlet boundary conditions by the least squares method
% Rectangle [0,c]x[0,b] or [-c,c]x[0,b]
% Boundary conditions: Dirichlet conditions

%clear all

% Number of points
%NX = 30;
%NY = 30;
% Spatial dimensions
b = 100;
c = 100;
a = 10; % length under the pressure p
DOMX = [-c,c];
DOMY = [0,b];
% DOMX = [0,c];
% DOMY = [0,b];
% Acting pressure
% p = 1e-01; pressure='Pa';
% p = 1e+06; pressure= 'MPa';
% Material constants
rat = 0.49;
E = 10e+09; 
% rat = 0.25;
% E = 200e+09; 
lambda = E*rat/((1+rat)*(1-2*rat));
mu = E/(2*(1+rat));
% Script name
ScriptName=join([mfilename,'_',pressure,'_',string(NX),string(NY),],'');

% Identity matrices
IdentityX = eye(NX);
IdentityY = eye(NY);

% Differentiation Matrices
Dx = diffmat(NX, 1,DOMX);
Dy = diffmat(NY, 1,DOMY);
Dx2 = diffmat(NX, 2,DOMX);
Dy2 = diffmat(NY, 2,DOMY);

DDy = kron(Dy,IdentityX);
DDx = kron(IdentityY,Dx);
DDy2 = kron(Dy2,IdentityX);
DDx2 = kron(IdentityY,Dx2);

% Boundary conditions
% Left
L=kron(IdentityY,IdentityX(1,:));
% Right
R = kron(IdentityY,IdentityX(end,:));
% Down
D = kron(IdentityY(end,:),IdentityX);
% Up
U = kron(IdentityY(1,:),IdentityX);

% Matrix

B = [zeros(size(U)),U,zeros(size(U))
    zeros(size(U)),zeros(size(U)),U
    zeros(size(D)),D,zeros(size(D))
    zeros(size(D)),zeros(size(D)),D
    L, zeros(size(L)),zeros(size(L))
    zeros(size(L(2:end-1,:))),zeros(size(L(2:end-1,:))),L(2:end-1,:)
    R, zeros(size(R)),zeros(size(R))
    zeros(size(R(2:end-1,:))),zeros(size(R(2:end-1,:))),R(2:end-1,:)]; 

Eq = [DDx,zeros(NX*NY),DDy
    zeros(NX*NY),DDy,DDx
    (1-rat).*DDy2-rat.*DDx2,(1-rat).*DDx2-rat.*DDy2, -DDx*DDy-DDy*DDx];
M=[B
   Eq];

% Right-hand side
% Analytical solution
Func_sigmax=@(x,y)  -(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
Func_sigmay=@(x,y)  -(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
Func_shearxy=@(x,y) -(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

sigmaxU = @(x)Func_sigmax(x,DOMY(1));
sigmayU = @(x)Func_sigmay(x,DOMY(1));
shearxyU = @(x)Func_shearxy(x,DOMY(1));

sigmaxL = @(y)Func_sigmax(DOMX(1),y);
sigmayL = @(y)Func_sigmay(DOMX(1),y);
shearxyL = @(y)Func_shearxy(DOMX(1),y);

sigmaxR = @(y)Func_sigmax(DOMX(2),y);
sigmayR = @(y)Func_sigmay(DOMX(2),y);
shearxyR = @(y)Func_shearxy(DOMX(2),y);

sigmaxD = @(x)Func_sigmax(x,DOMY(2));
sigmayD = @(x)Func_sigmay(x,DOMY(2));
shearxyD = @(x)Func_shearxy(x,DOMY(2));

sigxU = gridsample(sigmaxU,NX,DOMX);
sigyU = gridsample(sigmayU,NX,DOMX);
shearU = gridsample(shearxyU,NX,DOMX);

sigxD = gridsample(sigmaxD,NX,DOMX);
sigyD = gridsample(sigmayD,NX,DOMX);
shearD = gridsample(shearxyD,NX,DOMX);

sigxL = gridsample(sigmaxL,NY,DOMY);
sigyL = gridsample(sigmayL,NY,DOMY);
shearL = gridsample(shearxyL,NY,DOMY);

sigxR = gridsample(sigmaxR,NY,DOMY);
sigyR = gridsample(sigmayR,NY,DOMY);
shearR = gridsample(shearxyR,NY,DOMY);

rhs =[sigyU
    shearU
    sigyD
    shearD
    sigxL
    shearL(2:end-1,1)
    sigxR
    shearR(2:end-1,1)
    zeros(size(Eq,1),1)];

% Solution
tic;
sol = M\rhs;
time=toc;
fprintf('Elapsed time = %.3f seconds.\n', time)

% Reshape results

disx=sol(1:NX*NY);
MATsigx = transpose(reshape(disx,[NX,NY]));

disy=sol(NX*NY+1:2*NX*NY);
MATsigy = transpose(reshape(disy,[NX,NY]));

disxy=sol(2*NX*NY+1:3*NX*NY);
MATsh = transpose(reshape(disxy,[NX,NY]));

% Test + Plot results
%run('benchmarkS.m');

