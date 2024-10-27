% Domain boundary
function B = domain_boundary(NX,NY,domain)
    % Identity matrices
    IdentityX = eye(NX);
    IdentityY = eye(NY);
    
    % Boundary conditions
    % Left
    L=kron(IdentityY,IdentityX(1,:));
    % Right
    R = kron(IdentityY,IdentityX(end,:));
    % Down
    D = kron(IdentityY(end,:),IdentityX);
    % Up
    U = kron(IdentityY(1,:),IdentityX);
    
    zU=zeros(size(U));
    zR=zeros(size(R));
    zD=zeros(size(D));
    zL=zeros(size(L));
    
    % Matrix
    if domain=='R'
        B = [zU,U,zU
             zU,zU,U
             zD,D,zD
             zD,zD,D
             R,zR,zR
             zR(2:end-1,:),zR(2:end-1,:),R(2:end-1,:)];
    elseif domain=='L'
        B = [zU,U,zU
             zU,zU,U
             zD,D,zD
             zD,zD,D
             L, zL,zL
             zL(2:end-1,:),zL(2:end-1,:),L(2:end-1,:)];
    end
end