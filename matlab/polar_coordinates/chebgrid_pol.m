clear all

NR = 20;
NPhi = 20;

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

[r1,phi1]=meshgrid(r,phi);
[xx,yy] = pol2cart(phi1,r1);

figure('Name', 'Grid')
x=reshape(xx,[NR*NPhi,1]);
y=reshape(yy,[NR*NPhi,1]);
scatter(x,y,'k','.')
grid on
