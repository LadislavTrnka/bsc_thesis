% The script solving the Poisson equation with boundary conditions including discontinuities by the resampling technique from the Section 2.2.1
% Technique 1

NX = 40;NY = 40; % Number of points
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
% Up Sigma_y
xx = chebpts(NX);
jump = zeros(NX,1);
for i=1:NX
   if abs(xx(i))<0.2
      jump(i)=1;
   end
end

rhs=[jump
    zeros(NX+2*NY-4,1)
    zeros(size(P*DDx2,1),1)];

% Solution
tic;
sol = M\rhs;
time=toc;
fprintf('Elapsed time = %.3f seconds.\n', time)

% Reshape results and plot
MAT = transpose(reshape(sol,[NX,NY]));
[x, y] = chebpts2(NX, NY);
figure('Name', 'PoissonJump')
set(gcf,'units','centimeters','position',[0,0,20,20])
mesh(x,y,MAT,'edgecolor', 'k');
set(gca,'FontSize',10);
xlabel('$x$','interpreter','latex', 'FontWeight','bold','FontSize',12)
ylabel('$y$','interpreter','latex', 'FontWeight','bold','FontSize',12)
set(gca,'LooseInset',get(gca,'TightInset'));
exportgraphics(gcf,'..\Fig\PoissonJump.pdf','ContentType','image','Resolution',300)

