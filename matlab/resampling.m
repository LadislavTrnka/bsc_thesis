% To generate resampling matrices

function Res = resampling(NX,DOMX,NY,DOMY,n,m)

%Resampling X
x1 = chebpts(NX-n, DOMX,1);
x2 = chebpts(NX, DOMX,2);
resX = barymat(x1, x2);

%Resampling Y
y1 = chebpts(NY-m, DOMY,1);
y2 = chebpts(NY, DOMY,2);
resY = barymat(y1, y2);

R1 = kron(eye(NY-m),resX)*kron(resY,eye(NX));
R2 = kron(resY,eye(NX-n))*kron(eye(NY),resX);
Res = 1/2*(R1 + R2);
end

