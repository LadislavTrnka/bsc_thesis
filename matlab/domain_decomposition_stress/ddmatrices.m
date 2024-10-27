% Domain Decomposition matrices
function Eq = ddmatrices(NX,NY,DOMX,DOMY,letter)
    global rat
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
    
    Res=resampling(NX,DOMX,NY,DOMY,2,0);
    Res2=resampling(NX,DOMX,NY,DOMY,0,2);
    Res3=resampling(NX,DOMX,NY,DOMY,2,2);

%     Res=resampling(NX,DOMX,NY,DOMY,0,0);
%     Res2=resampling(NX,DOMX,NY,DOMY,0,0);
%     Res3=resampling(NX,DOMX,NY,DOMY,0,0);
    
    if letter=='A'
        P=blkdiag(Res2,Res3,Res);
    elseif letter=='B'
        P=blkdiag(Res3,Res,Res2);
    end
    
    % Matrix
    Eq = P*[DDx,zeros(NX*NY),DDy
        zeros(NX*NY),DDy,DDx
        (1-rat).*DDy2-rat.*DDx2,(1-rat).*DDx2-rat.*DDy2, -DDx*DDy-DDy*DDx];
end