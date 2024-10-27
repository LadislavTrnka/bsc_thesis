% The script solving the Poisson equation by the resampling technique from the Section 2.2.1
% Technique 1

% Number of points
NX = 20;NY = 20; 
IdentityY = eye(NY); IdentityX = eye(NX);

% Differentiation Matrices
DX2 = diffmat(NX, 2); DY2 = diffmat(NY, 2);
DDy2 = kron(DY2,IdentityX); DDx2 = kron(IdentityY,DX2);

% Resampling
x1 = chebpts(NX-2,1); y1 = chebpts(NY-2,1);
resX = barymat(x1, chebpts(NX,2)); % direction x
resY = barymat(y1, chebpts(NY,2)); % direction y
R1 = kron(eye(NY-2),resX)*kron(resY,IdentityX);
R2 = kron(resY,eye(NX-2))*kron(IdentityY,resX);
P = 1/2*(R1 + R2);

% Boundary conditions
L = kron(IdentityY,IdentityX(1,:)); % Left
R = kron(IdentityY,IdentityX(end,:)); % Right
D = kron(IdentityY(end,:),IdentityX); % Down
U = kron(IdentityY(1,:),IdentityX); % Up

% Matrix
B = [U
   D
   L(2:end-1,:)
   R(2:end-1,:)]; 
M = [B
   P*(DDx2+DDy2)];

% Right-hand side
Func_F = @(x,y) -8*pi^2*sin(2*pi*x)*sin(2*pi*y); F = Func_F(x1,y1');
rhs=[zeros(2*NX+2*NY-4,1)
    reshape(F,[],1)];

% Solution
sol = M\rhs;

% Reshape results
MAT = transpose(reshape(sol,[NX,NY]));
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
