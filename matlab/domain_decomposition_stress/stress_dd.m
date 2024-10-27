% Stress_DomainDecomposition - resampling

clear all

% Number of points in domains
NX1 = 50;
NX2 = 50;
NY = 50;
% Spatial dimensions
b = 100;
c = 100;
global a
a = 10; % length under the pressure p
DOMY = [5,10];
DOMX1 = [-8,0];
DOMX2 = [0,8];
% Acting pressure
global p
p = 1e-01; pressure='Pa';
% p = 1e+06; pressure= 'MPa';

% Material constants
global rat
rat = 0.49;

% Equations and boundary conditions
M1 = ddmatrices(NX1,NY,DOMX1,DOMY,'A');
M2 = ddmatrices(NX2,NY,DOMX2,DOMY,'B');
B1 = domain_boundary(NX1,NY,'L');
B2 = domain_boundary(NX2,NY,'R');
SB = shared_boundary(NX1,DOMX1,NX2,DOMX2,NY);

% Matrix
M = [B1, zeros(size(B1,1),size(B2,2))
    zeros(size(B2,1),size(B1,2)),B2
    SB
    M1, zeros(size(M1,1),size(M2,2))
    zeros(size(M2,1),size(M1,2)),M2]; 

% Right hand side
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

MATsigx1 = transpose(reshape(sol(1:NX1*NY),[NX1,NY]));
MATsigy1 = transpose(reshape(sol(NX1*NY+1:2*NX1*NY),[NX1,NY]));
MATsh1 = transpose(reshape(sol(2*NX1*NY+1:3*NX1*NY),[NX1,NY]));

MATsigx2 = transpose(reshape(sol(3*NX1*NY+1:3*NX1*NY+NX2*NY),[NX2,NY]));
MATsigy2 = transpose(reshape(sol(3*NX1*NY+NX2*NY+1:3*NX1*NY+2*NX2*NY),[NX2,NY]));
MATsh2 = transpose(reshape(sol(3*NX1*NY+2*NX2*NY+1:3*NX1*NY+3*NX2*NY),[NX2,NY]));

% % Test + Plot results
% run('benchmarkS.m');
