% To compute stress formulation with Dirichlet boundary conditions by the moving technique
% Rectangle [0,c]x[0,b] or [-c,c]x[0,b]
% Boundary conditions: Dirichlet conditions

%clear all

% Number of points
%NX = 40;
%NY = 40;
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
% rat = 0.25;

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
% we have 3 options:
%letter='A';
%letter='B';
if letter=='A'
    % A - oscillations in shear
Eq = {zeros(NX*NY),DDy,DDx
    (1-rat).*DDy2-rat.*DDx2,(1-rat).*DDx2-rat.*DDy2, -2.*DDx*DDy
    DDx,zeros(NX*NY),DDy};
elseif letter=='B'
    % B - oscillations in sigma_x and sigma_y
Eq = {(1-rat).*DDy2-rat.*DDx2,(1-rat).*DDx2-rat.*DDy2, -2.*DDx*DDy
     DDx,zeros(NX*NY),DDy
     zeros(NX*NY),DDy,DDx};
end
% Large condition number e+18 for 30x30
% Eq = {zeros(NX*NY),DDy,DDx
%     DDx,zeros(NX*NY),DDy
%     (1-rat).*DDy2-rat.*DDx2,(1-rat).*DDx2-rat.*DDy2, -2.*DDx*DDy};

% Script name
ScriptName=join([mfilename,'_',pressure,'_',string(NX),string(NY),letter],'');

for i=1:size(Eq,1)
    for j=1:size(Eq,2)
        Eq{i,j}=removing_stress(cell2mat(Eq(i,j)),NX,NY,1,i);
    end
end
    
C = cell2mat(Eq);

% Right-hand side
% Analytical solution
Func_sigmax=@(x,y)  -(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
Func_sigmay=@(x,y)  -(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
Func_shearxy=@(x,y) -(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

sigmayU = @(x)Func_sigmay(x,DOMY(1));
shearxyU = @(x)Func_shearxy(x,DOMY(1));

sigmaxL = @(y)Func_sigmax(DOMX(1),y);
shearxyL = @(y)Func_shearxy(DOMX(1),y);

sigmaxR = @(y)Func_sigmax(DOMX(2),y);
shearxyR = @(y)Func_shearxy(DOMX(2),y);

sigmayD = @(x)Func_sigmay(x,DOMY(2));
shearxyD = @(x)Func_shearxy(x,DOMY(2));

sigyU = gridsample(sigmayU,NX,DOMX);
shearU = gridsample(shearxyU,NX,DOMX);

sigyD = gridsample(sigmayD,NX,DOMX);
shearD = gridsample(shearxyD,NX,DOMX);

sigxL = gridsample(sigmaxL,NY,DOMY);
shearL = gridsample(shearxyL,NY,DOMY);

sigxR = gridsample(sigmaxR,NY,DOMY);
shearR = gridsample(shearxyR,NY,DOMY);

% Organize the boundary vector sigmx_B, sigmy_B and shear_B

P = per(NY,NX);

sigmx_B = P*[sigxL; zeros(NY*(NX-2),1); sigxR];

sigmy_B = zeros(NY*NX,1);
sigmy_B(1:NX,1)=sigyU;
sigmy_B(end-NX+1:end,1)=sigyD;

shear_B = P*[shearL; zeros(NY*(NX-2),1); shearR];
shear_B(1:NX,1)=shearU;
shear_B(end-NX+1:end,1)=shearD;

rhs=-C*[sigmx_B;sigmy_B;shear_B];

% Matrix M
for i=1:size(Eq,1)
    for j=1:size(Eq,2)
        Eq{i,j}=removing_stress(cell2mat(Eq(i,j)),NX,NY,0,j);
    end
end
    
M = cell2mat(Eq);

% Solution
tic;
sol = M\rhs;
time=toc;
fprintf('Elapsed time = %.3f seconds.\n', time)

% Reshape results

disx=sol(1:(NX-2)*NY);
MATsigmax = transpose(reshape(disx,[NX-2,NY]));

MATsigx = zeros(NY,NX);
MATsigx(1:end,1)=sigxL;
MATsigx(1:end,end)=sigxR;
MATsigx(:, 2 : (end - 1))=MATsigmax;

disy=sol((NX-2)*NY+1:(NX-2)*NY+NX*(NY-2));
MATsigmay = transpose(reshape(disy,[NX,NY-2]));

MATsigy = zeros(NY,NX);
MATsigy(1,1:end)=sigyU;
MATsigy(end,1:end)=sigyD;
MATsigy(2 : (end - 1), :)=MATsigmay;

disxy=sol((NX-2)*NY+NX*(NY-2)+1:end);
MATshear = transpose(reshape(disxy,[NX-2,NY-2]));

MATsh = zeros(NY,NX);
MATsh(1:end,1)=shearL;
MATsh(1:end,end)=shearR;
MATsh(1,1:end)=shearU;
MATsh(end,1:end)=shearD;
MATsh(2 : (end - 1), 2 : (end - 1))=MATshear;

% Test + Plot results
%run('benchmarkS.m');
