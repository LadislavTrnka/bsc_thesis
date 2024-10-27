% Domain right hand side
function rhs = domain_rhs(NX,NY,DOMX,DOMY,Boundary)
    global p a
    % Right hand side
    % Analytical solution
    Func_sigmax=@(x,y)  -(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
    Func_sigmay=@(x,y)  -(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
    Func_shearxy=@(x,y) -(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

    %sigmaxU = @(x)Func_sigmax(x,DOMY(1));
    sigmayU = @(x)Func_sigmay(x,DOMY(1));
    shearxyU = @(x)Func_shearxy(x,DOMY(1));

    sigmaxL = @(y)Func_sigmax(DOMX(1),y);
    %sigmayL = @(y)Func_sigmay(DOMX(1),y);
    shearxyL = @(y)Func_shearxy(DOMX(1),y);

    sigmaxR = @(y)Func_sigmax(DOMX(2),y);
    %sigmayR = @(y)Func_sigmay(DOMX(2),y);
    shearxyR = @(y)Func_shearxy(DOMX(2),y);

    %sigmaxD = @(x)Func_sigmax(x,DOMY(2));
    sigmayD = @(x)Func_sigmay(x,DOMY(2));
    shearxyD = @(x)Func_shearxy(x,DOMY(2));

    %sigxU = gridsample(sigmaxU,NX,DOMX);
    sigyU = gridsample(sigmayU,NX,DOMX);
    shearU = gridsample(shearxyU,NX,DOMX);

    %sigxD = gridsample(sigmaxD,NX,DOMX);
    sigyD = gridsample(sigmayD,NX,DOMX);
    shearD = gridsample(shearxyD,NX,DOMX);

    sigxL = gridsample(sigmaxL,NY,DOMY);
    %sigyL = gridsample(sigmayL,NY,DOMY);
    shearL = gridsample(shearxyL,NY,DOMY);

    sigxR = gridsample(sigmaxR,NY,DOMY);
    %sigyR = gridsample(sigmayR,NY,DOMY);
    shearR = gridsample(shearxyR,NY,DOMY);
    
    if Boundary=='L'   
        rhs =[sigyU
            shearU
            sigyD
            shearD
            sigxL
            shearL(2:end-1,1)];
        rhs(isnan(rhs))=-p;
    elseif Boundary=='R'   
        rhs =[sigyU
            shearU
            sigyD
            shearD
            sigxR
            shearR(2:end-1,1)];
        rhs(isnan(rhs))=0;
    end
end