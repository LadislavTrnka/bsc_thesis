% To compute displacement formulation with Dirichlet boundary conditions by the resampling technique
% Rectangle [0,c]x[0,b]
% Boundary conditions: Dirichlet conditions

%clear all

% Number of points
%NX = 30;
%NY = 30;
% Spatial dimensions
b = 100;
c = 100;
a = 10; % length under the pressure p
DOMX = [0,c];
DOMY = [0,b];
% Acting pressure
%p = 1e-01; pressure='Pa';
%p = 1e+06; pressure= 'MPa';
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
L(1,:)=[];
L(end,:)=[];
% Right
R = kron(IdentityY,IdentityX(end,:));
R(1,:)=[];
R(end,:)=[];
% Down
D = kron(IdentityY(end,:),IdentityX);
% Up
U = kron(IdentityY(1,:),IdentityX);

% Resampling
Res=resampling(NX,DOMX,NY,DOMY,2,2);
P=blkdiag(Res,Res);

% Matrix

B = [U,zeros(size(U))
    zeros(size(U)),U
    L,zeros(size(L))
    zeros(size(L)),L
    R,zeros(size(R))
    zeros(size(R)),R
    D,zeros(size(D))
    zeros(size(D)),D]; 

% Due to large dimension of matherial constants, we add factor 1e-9 
scale=1e-9;
Eq = scale.*[(lambda+2*mu).*DDx2+mu*DDy2,(lambda+mu).*1/2*(DDx*DDy+DDy*DDx)
    (lambda+mu).*1/2*(DDx*DDy+DDy*DDx), (lambda+2*mu).*DDy2+mu*DDx2];
M=[B
   P*Eq];

% Right-hand side
% Analytical solution
Func_u=@(x,y)(1+rat)./(pi.*E)*p.*(-(2.*rat-1).*(a-x).*atan((a-x)./y)+(2.*rat-1).*(a+x).*atan((a+x)./y)+(rat-1).*y.*log((y.^2+(a-x).^2)./(y.^2+(a+x).^2)));
zero_level=@(x,y)  (2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2));
Func_v=@(x,y) (1+rat)./(pi.*E)*p.*((2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2))-zero_level(c,b));

FuU = @(x)Func_u(x,DOMY(1));
FvU = @(x)Func_v(x,DOMY(1));

FuL = @(y)Func_u(DOMX(1),y);
FvL = @(y)Func_v(DOMX(1),y);

FuR = @(y)Func_u(DOMX(2),y);
FvR = @(y)Func_v(DOMX(2),y);

FuD = @(x)Func_u(x,DOMY(2));
FvD = @(x)Func_v(x,DOMY(2));

uU = gridsample(FuU,NX,DOMX);
vU = gridsample(FvU,NX,DOMX);

uD = gridsample(FuD,NX,DOMX);
vD = gridsample(FvD,NX,DOMX);

uR = gridsample(FuR,NY,DOMY);
vR = gridsample(FvR,NY,DOMY);

uL = gridsample(FuL,NY,DOMY);
vL = gridsample(FvL,NY,DOMY);

% Right hand side

rhs = [uU
    vU
    uL(2:end-1,1)
    vL(2:end-1,1)
    uR(2:end-1,1)
    vR(2:end-1,1)
    uD
    vD
    zeros(size(P*Eq,1),1)];

% Solution
tic;
sol = M\rhs;
time=toc;
fprintf('Elapsed time = %.3f seconds.\n', time)

% Reshape results

disu=sol(1:NX*NY);
MATu =transpose(reshape(disu,[NX,NY]));

disv=sol(NX*NY+1:end);
MATv =transpose(reshape(disv,[NX,NY]));

% Test + Plot results
%run('benchmarkD.m');