% The script solving the Poisson equation by the moving technique from the Section 2.2.2
% Technique 2

% Number of points
NX = 20; NY = 20; 

% Differentiation Matrices
DX2 = diffmat(NX, 2); DY2 = diffmat(NY, 2);
DDy2 = kron(DY2,eye(NX)); DDx2 = kron(eye(NY),DX2);

% Position of boundary conditions
U = 1:NX; % up 
R = NX:NX:NX*NY; % right
L = 1:NX:NX*NY; % left
D = NX*(NY-1)+1:NX*NY; % down
k = [L R U D]; % all positions

% Matrix C
C = DDx2+DDy2; C(k,:)=[];

% Right-hand side
Func_F = @(x,y) -8*pi^2*sin(2*pi*x)*sin(2*pi*y);
F = Func_F(chebpts(NX,2),chebpts(NY,2)'); F = reshape(F,[],1);
F(k,:) = []; rhs = F; %u_b=0

% Matrix M
C(:,k)=[]; M = C;

% Solution
sol = M\rhs;

% Reshape results
MAT=zeros(NX,NY);
MAT(2:end-1,2:end-1)=transpose(reshape(sol,[NX-2,NY-2]));
u = chebfun2(MAT);
plot(u)

[x, y] = chebpts2(NX, NY);
z= sin(2.*pi.*x).*sin(2.*pi.*y);
error = abs(z-MAT);
Name = ["u"];
NormInf = [norm(error,'Inf')];
MaxError = [max(max(abs(error)))];
table(NormInf,MaxError, 'RowNames',Name)
figure('Name','error')
mesh(x,y,error);
