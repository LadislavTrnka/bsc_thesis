% Displacement_DomainDecomposition (2) - resampling

clear all

% Number of points in domains
NX1 = 80;
NX2 = 80;
NY = 80;
% Spatial dimensions
b = 100;
c = 100;
global a p rat E lambda mu scale
a = 10; % length under the pressure p
DOMY = [0,b];
DOMX1 = [0,a];
DOMX2 = [a,c];
% Acting pressure = size of the jump - larger pressure, bigger oscillations in solution
p = 1e-01; pressure='Pa';
% p = 1e+06; pressure= 'MPa';

% Material constants
rat = 0.49;
E = 10e+09; 
lambda = E*rat/((1+rat)*(1-2*rat));
mu = E/(2*(1+rat));


% Equations and boundary conditions
scale=1e-9;
M1 = ddmatrices(NX1,NY,DOMX1,DOMY);
M2 = ddmatrices(NX2,NY,DOMX2,DOMY);
B1 = domain_boundary(NX1,NY,DOMX1,DOMY,'L');
B2 = domain_boundary(NX2,NY,DOMX2,DOMY,'R');
SB = shared_boundary(NX1,DOMX1,NX2,DOMX2,NY); 

% Matrix
M = [B1, zeros(size(B1,1),size(B2,2))
    zeros(size(B2,1),size(B1,2)),B2
    SB
    M1, zeros(size(M1,1),size(M2,2))
    zeros(size(M2,1),size(M1,2)),M2]; 

% Right hand side
global zl
zero_level=@(x,y)  (2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2));
zl=zero_level(0,b);

rhs1=domain_rhs(NX1,NY,DOMX1,DOMY,'L');
rhs2=domain_rhs(NX2,NY,DOMX2,DOMY,'R');
rhs = [rhs1
    rhs2
    zeros(size(SB,1),1)
    zeros(size(M1,1),1)
    zeros(size(M2,1),1)];

% Solution
tic;
sol = M\rhs;
time=toc;
fprintf('Elapsed time = %.3f seconds.\n', time)

% Reshape results

MATu1 =transpose(reshape(sol(1:NX1*NY),[NX1,NY]));
MATv1 =transpose(reshape(sol(NX1*NY+1:2*NX1*NY),[NX1,NY]));

MATu2 =transpose(reshape(sol(2*NX1*NY+1:2*NX1*NY+NX2*NY),[NX2,NY]));
MATv2 =transpose(reshape(sol(2*NX1*NY+NX2*NY+1:2*NX1*NY+2*NX2*NY),[NX2,NY]));

% % Test + Plot results
% run('benchmarkS.m');
