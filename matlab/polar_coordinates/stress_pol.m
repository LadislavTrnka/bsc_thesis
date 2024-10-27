clear all

NR = 90;
NPhi = 90;

a = 10; % length under the pressure p
c = 100;
%a = 100;
%c = 10;
b = c;
p = 1;

DOMR = [0,c];
DOMPhi = [0,pi];
r=chebpts(NR, DOMR, 1);
phi=chebpts(NPhi, DOMPhi);

FUN_sigmax=@(x,y)  -(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
FUN_sigmay=@(x,y)  -(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
FUN_shearxy=@(x,y) -(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

[r1,phi1]=meshgrid(r,phi);
[x,y] = pol2cart(phi1,r1);

mesh(x,y,FUN_sigmax(x,y))
mesh(x,y,FUN_sigmay(x,y))
mesh(x,y,FUN_shearxy(x,y))

sigR=FUN_sigmax(x,y).*(cos(phi1).^2)+FUN_sigmay(x,y).*(sin(phi1).^2)+2.*FUN_shearxy(x,y).*(sin(phi1).*cos(phi1));
mesh(x,y,sigR)
sigPhi=FUN_sigmax(x,y).*(sin(phi1).^2)+FUN_sigmay(x,y).*(cos(phi1).^2)-2.*FUN_shearxy(x,y).*(sin(phi1).*cos(phi1));
mesh(x,y,sigPhi)
shear=(FUN_sigmay(x,y)-FUN_sigmax(x,y)).*(sin(phi1).*cos(phi1))+FUN_shearxy(x,y).*((cos(phi1).^2)-(sin(phi1).^2));
mesh(x,y,shear)

scatter(reshape(x,[NR*NPhi,1]),reshape(y,[NR*NPhi,1]),'filled')
