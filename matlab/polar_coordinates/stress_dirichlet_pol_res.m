clear all

NR = 30;
NPhi = 30;

a = 10; % length under the pressure p
c = 100;
%a = 100;
%c = 10;
b = c;
p = 1;

DOMR = [0,c];
DOMPhi = [0,pi];

rat = 0.49;

DR = diffmat(NR, 1,DOMR,'chebkind1');
DPhi = diffmat(NPhi, 1,DOMPhi);

DR2 = diffmat(NR, 2,DOMR,'chebkind1');
DPhi2 = diffmat(NPhi, 2,DOMPhi);

DDR = kron(eye(NPhi),DR);
DDPhi = kron(DPhi,eye(NR));

DDR2 = kron(eye(NPhi),DR2);
DDPhi2 = kron(DPhi2,eye(NR));


%Hranice
IdentityR = eye(NR);
IdentityPhi = eye(NPhi);

U1= kron(IdentityPhi(1,:),IdentityR);
U2= kron(IdentityPhi(end,:),IdentityR);

RadiusOUT=kron(IdentityPhi,IdentityR(end,:));
RadiusShOUT=RadiusOUT(2:end-1,:);

RadiusIN=kron(IdentityPhi,IdentityR(1,:));
RadiusShIN=RadiusIN(2:end-1,:);
%Matrix
r=chebpts(NR, DOMR, 1);
phi=chebpts(NPhi, DOMPhi);
rr=kron(IdentityPhi,diag(1./r));

% A=DDR2+rr*DDR+rr.^2*DDPhi2;
A=(1-rat).*rr.^2*DDPhi2-rat.*DDR2-(1+rat).*rr*DDR;
B=-rat.*rr.^2*DDPhi2+(1-rat).*DDR2+(2-rat).*rr*DDR;
C=-rr*(DDPhi*DDR+DDR*DDPhi)-2.*rr.^2*DDPhi;

R2=resampling_pol(NR,DOMR,NPhi,DOMPhi,2,0);
R3=resampling_pol(NR,DOMR,NPhi,DOMPhi,0,2);
R1=resampling_pol(NR,DOMR,NPhi,DOMPhi,2,2);

M = [zeros(size(U1)),U1,zeros(size(U1))
    zeros(size(U1)),zeros(size(U1)),U1
    zeros(size(U2)),U2,zeros(size(U2))
    zeros(size(U2)),zeros(size(U2)),U2
    RadiusOUT, zeros(size(RadiusOUT)),zeros(size(RadiusOUT))
    zeros(size(RadiusShOUT)),zeros(size(RadiusShOUT)),RadiusShOUT
    RadiusIN, zeros(size(RadiusIN)),zeros(size(RadiusIN))
    zeros(size(RadiusShIN)),zeros(size(RadiusShIN)),RadiusShIN
    R1*(DDR+rr),-R1*rr,R1*(rr*DDPhi)
    zeros(size(R2*DDPhi)),R2*(rr*DDPhi), R2*(DDR+2.*rr)
    R3*A,R3*B, R3*C];

%ConditionNumber=cond(M);
rank(M)
%Right hand side
RHSsigPhi = zeros(NR,1);
for i=1:NR
   if r(i)<=a
        RHSsigPhi(i)=-p;
   end
end

%Analytical solution
[r0,phi0]=meshgrid(r(end),phi);
[x0,y0] = pol2cart(phi0,r0);

FUN_sigmax=@(x,y)  -(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
FUN_sigmay=@(x,y)  -(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
FUN_shearxy=@(x,y) -(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

RHSsigR=FUN_sigmax(x0,y0).*(cos(phi0).^2)+FUN_sigmay(x0,y0).*(sin(phi0).^2)+2.*FUN_shearxy(x0,y0).*(sin(phi0).*cos(phi0));
RHSshear=(FUN_sigmay(x0,y0)-FUN_sigmax(x0,y0)).*(sin(phi0).*cos(phi0))+FUN_shearxy(x0,y0).*((cos(phi0).^2)-(sin(phi0).^2));

rhs =[RHSsigPhi
    zeros(size(U1,1),1)
    RHSsigPhi
    zeros(size(U2,1),1)
    RHSsigR
    RHSshear(2:end-1)
    -ones(size(RadiusIN,1),1)
    zeros(size(RadiusShIN,1),1)
    zeros(size(R1*DDPhi,1),1)
    zeros(size(R2*DDPhi,1),1)
    zeros(size(R3*A,1),1)];

%System of linear equations
tic
sol = M\rhs;
toc

disx=sol(1:NR*NPhi);
MATsigR = transpose(reshape(disx,[NR,NPhi]));

disy=sol(NR*NPhi+1:2*NR*NPhi);
MATsigPhi = transpose(reshape(disy,[NR,NPhi]));

disxy=sol(2*NR*NPhi+1:3*NR*NPhi);
MATshRPhi = transpose(reshape(disxy,[NR,NPhi]));

[Mr,Mphi]=meshgrid(r,phi);
[x,y] = pol2cart(Mphi,Mr);

figure(1)
subplot(1,3,1)
mesh(x,y,MATsigR)
subplot(1,3,2)
mesh(x,y,MATsigPhi)
subplot(1,3,3)
mesh(x,y,MATshRPhi)

MATsigx=MATsigR.*(cos(Mphi).^2)+MATsigPhi.*(sin(Mphi).^2)-2.*MATshRPhi.*(sin(Mphi).*cos(Mphi));
MATsigy=MATsigR.*(sin(Mphi).^2)+MATsigPhi.*(cos(Mphi).^2)+2.*MATshRPhi.*(sin(Mphi).*cos(Mphi));
MATsh=(MATsigR-MATsigPhi).*(sin(Mphi).*cos(Mphi))+MATshRPhi.*((cos(Mphi).^2)-(sin(Mphi).^2));

