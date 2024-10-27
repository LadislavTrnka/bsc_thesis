% Common boundary
function B = shared_boundary(NXL,DOMXL,NXR,DOMXR,NY)
    % Identity matrices
    IdentityXL = eye(NXL);
    IdentityXR = eye(NXR);
    IdentityY = eye(NY);
    % Differentiation Matrices
    DxL = diffmat(NXL, 1,DOMXL);
    DxR = diffmat(NXR, 1,DOMXR);

    % Boundary conditions
    % Left
    L = kron(IdentityY,IdentityXR(1,:));
    % Right
    R = kron(IdentityY,IdentityXL(end,:));
    % Left - First derivative
    DL = kron(IdentityY,DxR(1,:));
    % Right - First derivative
    DR = kron(IdentityY,DxL(end,:));
    
    zR=zeros(size(R));
    zL=zeros(size(L));
    
    % Matrix
    BL = [L(2:end-1,:),zL(2:end-1,:)
        zL(2:end-1,:),L(2:end-1,:)]; 
    BR = [R(2:end-1,:),zR(2:end-1,:)
          zR(2:end-1,:),R(2:end-1,:)]; 
    BDL = [DL(2:end-1,:),zL(2:end-1,:)
        zL(2:end-1,:),DL(2:end-1,:)]; 
    BDR = [DR(2:end-1,:),zR(2:end-1,:)
          zR(2:end-1,:),DR(2:end-1,:)];     
    B = [BR,-BL
        BDR,-BDL];
end