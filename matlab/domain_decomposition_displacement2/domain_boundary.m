% Domain boundary
function B = domain_boundary(NX,NY,DOMX,DOMY,domain)
    global lambda mu scale
    % Identity matrices
    IdentityX = eye(NX);
    IdentityY = eye(NY);
    % Differentiation Matrices
    Dx = diffmat(NX, 1,DOMX);
    Dy = diffmat(NY, 1,DOMY);
    
    % Up Sigma_y
    USy_u = scale.*(lambda.*kron(IdentityY(1,:),diffmat([NX,NX], 1,DOMX)));
    USy_v = scale.*((lambda+2*mu).*kron(Dy(1,:),IdentityX));

    % Up Shear
    UShear_u = kron(Dy(1,:),IdentityX);
    UShear_v = kron(IdentityY(1,:),diffmat([NX,NX], 1,DOMX));    

    % Down
    D = kron(IdentityY(end,:),IdentityX);    

    % Domain
    if domain=='L'
        % Left
        LDx = kron(IdentityY,Dx(1,:));
        LDx(1,:)=[];
        LDx(end,:)=[];
        L=kron(IdentityY,IdentityX(1,:));
        L(1,:)=[];
        L(end,:)=[];

        B = [USy_u,USy_v
            UShear_u,UShear_v
            zeros(size(LDx)),LDx
            L,zeros(size(L))
            D,zeros(size(D))
            zeros(size(D)),D];
    elseif domain=='R'
        % Right
        R = kron(IdentityY,IdentityX(end,:));
        R(1,:)=[];
        R(end,:)=[];

        B = [USy_u,USy_v
            UShear_u,UShear_v
            R,zeros(size(R))
            zeros(size(R)),R
            D,zeros(size(D))
            zeros(size(D)),D];
    end
end