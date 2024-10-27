function Res = resampling_pol(NR,DOMR,NPhi,DOMPhi,n,m)

%Resampling X
[r, ~, rBW]=chebpts(NR, DOMR, 1);
r2 = chebpts(NR-n, DOMR,1);
resX = barymat(r2, r, rBW);

%RESAMPLINX Y
phi=chebpts(NPhi, DOMPhi);
phi2 = chebpts(NPhi-m, DOMPhi,1);
resY = barymat(phi2, phi);

R1 = kron(eye(NPhi-m),resX)*kron(resY,eye(NR));
R2 = kron(resY,eye(NR-n))*kron(eye(NPhi),resX);
Res = 1/2*(R1 + R2);
end
