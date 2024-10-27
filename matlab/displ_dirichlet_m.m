% To compute displacement formulation with Dirichlet boundary conditions by the moving technique
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

% Matrix C
Eq = {(lambda+2*mu).*DDx2+mu*DDy2,(lambda+mu).*1/2*(DDx*DDy+DDy*DDx)
    (lambda+mu).*1/2*(DDx*DDy+DDy*DDx), (lambda+2*mu).*DDy2+mu*DDx2};

% Removing rows of Equations
for i=1:size(Eq,1)
    for j=1:size(Eq,2)
        Eq{i,j}=removing(cell2mat(Eq(i,j)),NX,NY,1);
    end
end
C = cell2mat(Eq);

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

% Organize the boundary vector u_B and u_B
P = per(NY,NX);
u_B = P*[uL; zeros(NY*(NX-2),1); uR];
u_B(1:NX,1)=uU;
u_B(end-NX+1:end,1)=uD;

v_B = P*[vL; zeros(NY*(NX-2),1); vR];
v_B(1:NX,1)=vU;
v_B(end-NX+1:end,1)=vD;

rhs=-C*[u_B;v_B];

% Matrix M
for i=1:size(Eq,1)
    for j=1:size(Eq,2)
        Eq{i,j}=removing(cell2mat(Eq(i,j)),NX,NY,0);
    end
end
M = cell2mat(Eq);

% Solution
tic;
sol = M\rhs;
time=toc;
fprintf('Elapsed time = %.3f seconds.\n', time)

% Reshape results

disu=sol(1:(NX-2)*(NY-2));
MATu = zeros(NY,NX);
MATu(1:end,1)=uL;
MATu(1:end,end)=uR;
MATu(end,1:end)=uD;
MATu(1,1:end)=uU;
MATu(2 : (end - 1), 2 : (end - 1))=transpose(reshape(disu,[NX-2,NY-2]));

disv=sol((NX-2)*(NY-2)+1:end);
MATv = zeros(NY,NX);
MATv(1:end,1)=vL;
MATv(1:end,end)=vR;
MATv(1,1:end)=vU;
MATv(end,1:end)=vD;
MATv(2 : (end - 1), 2 : (end - 1))=transpose(reshape(disv,[NX-2,NY-2]));

% Test + Plot results
%run('benchmarkD.m');