% Domain right hand side

function rhs = DomainRHS(NX,NY,DOMX,DOMY,Boundary)
    global p a scale rat E zl
    % Up Sigma_y
    xx = chebpts(NX, DOMX);
    sigyU = zeros(NX,1);
    for i=1:NX
       if abs(xx(i))<=a
            sigyU(i)=-p;
       end
    end

    % Analytical solution
    Func_u=@(x,y)(1+rat)./(pi.*E)*p.*(-(2.*rat-1).*(a-x).*atan((a-x)./y)+(2.*rat-1).*(a+x).*atan((a+x)./y)+(rat-1).*y.*log((y.^2+(a-x).^2)./(y.^2+(a+x).^2)));
    Func_v=@(x,y) (1+rat)./(pi.*E)*p.*((2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2))-zl);

%     FuU = @(x)Func_u(x,DOMY(1));
%     FvU = @(x)Func_v(x,DOMY(1));

    FuL = @(y)Func_u(DOMX(1),y);
    FvL = @(y)Func_v(DOMX(1),y);

    FuR = @(y)Func_u(DOMX(2),y);
    FvR = @(y)Func_v(DOMX(2),y);

    FuD = @(x)Func_u(x,DOMY(2));
    FvD = @(x)Func_v(x,DOMY(2));

%     uU = gridsample(FuU,NX,DOMX);
%     vU = gridsample(FvU,NX,DOMX);

    uD = gridsample(FuD,NX,DOMX);
    vD = gridsample(FvD,NX,DOMX);

    uR = gridsample(FuR,NY,DOMY);
    vR = gridsample(FvR,NY,DOMY);

    uL = gridsample(FuL,NY,DOMY);
    vL = gridsample(FvL,NY,DOMY);

    if Boundary=='R'   
        rhs = [scale.*sigyU
            zeros(NX,1)
            zeros(NY-2,1)
            zeros(NY-2,1)
            uR(2:end-1,1)
            vR(2:end-1,1)
            uD
            vD];
    elseif Boundary=='L'   
        rhs = [scale.*sigyU
            zeros(NX,1)
            uL(2:end-1,1)
            vL(2:end-1,1)
            uD
            vD];
    end
end