% Common boundary

function B = SharedBoundary(NXL,NXR,NY)
    % Identity matrices
    IdentityXL = eye(NXL);
    IdentityXR = eye(NXR);
    IdentityY = eye(NY);
    
    % Boundary conditions
    % Left
    L = kron(IdentityY,IdentityXR(1,:));
    % Right
    R = kron(IdentityY,IdentityXL(end,:));
    
    zR=zeros(size(R));
    zL=zeros(size(L));
    
    % Matrix
    BL = [L(2:end-1,:),zL(2:end-1,:)
        zL(2:end-1,:),L(2:end-1,:)]; 
    BR = [-R(2:end-1,:),zR(2:end-1,:)
          zR(2:end-1,:),-R(2:end-1,:)]; 
    B = [BL,BR];
end