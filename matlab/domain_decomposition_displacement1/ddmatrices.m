% Domain Decomposition matrices

function Eq = DomainDecompositionMatrices(NX,NY,DOMX,DOMY,domain)
    global lambda mu
    scale=1e-9;
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
    
    % Resampling
    if domain=='L'
        Res=resampling(NX,DOMX,NY,DOMY,2,2);
        P=blkdiag(Res,Res);
    elseif domain=='R'
        Res=resampling(NX,DOMX,NY,DOMY,2,2);
        P=blkdiag(Res,Res);
    end
    % Matrix
    Eq = scale.*P*[(lambda+2*mu).*DDx2+mu*DDy2,(lambda+mu).*1/2*(DDx*DDy+DDy*DDx)
        (lambda+mu).*1/2*(DDx*DDy+DDy*DDx), (lambda+2*mu).*DDy2+mu*DDx2];
end