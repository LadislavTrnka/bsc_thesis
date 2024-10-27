%Tests + plots: displacement (1)

[x1, y1] = chebpts2(NX1, NY, [DOMX1(1) DOMX1(2) DOMY(1) DOMY(2)]);
[x2, y2] = chebpts2(NX2, NY, [DOMX2(1) DOMX2(2) DOMY(1) DOMY(2)]);

zu1=(1+rat)./(pi.*E)*p.*(-(2.*rat-1).*(a-x1).*atan((a-x1)./y1)+(2.*rat-1).*(a+x1).*atan((a+x1)./y1)+(rat-1).*y1.*log((y1.^2+(a-x1).^2)./(y1.^2+(a+x1).^2)));
zv1=(1+rat)./(pi.*E)*p.*((2.*rat-1).*y1.*(atan((a-x1)./y1)+atan((a+x1)./y1))+(rat-1).*((a-x1).*log(y1.^2+(a-x1).^2)+(a+x1).*log(y1.^2+(a+x1).^2))-zl);

zu2=(1+rat)./(pi.*E)*p.*(-(2.*rat-1).*(a-x2).*atan((a-x2)./y2)+(2.*rat-1).*(a+x2).*atan((a+x2)./y2)+(rat-1).*y2.*log((y2.^2+(a-x2).^2)./(y2.^2+(a+x2).^2)));
zv2=(1+rat)./(pi.*E)*p.*((2.*rat-1).*y2.*(atan((a-x2)./y2)+atan((a+x2)./y2))+(rat-1).*((a-x2).*log(y2.^2+(a-x2).^2)+(a+x2).*log(y2.^2+(a+x2).^2))-zl);

x=[x2,x1];
y=[y2,y1];

MATu=[MATu2,MATu1];
MATv=[MATv2,MATv1];
zu=[zu2,zu1];
zv=[zv2,zv1];

%Error matrices
erroru = zu-MATu;
errorv = zv-MATv;

subplot(2,2,1);
mesh(x,y,MATu)
subplot(2,2,2);
mesh(x,y,MATv)
subplot(2,2,3);
mesh(x,y,erroru)
subplot(2,2,4);
mesh(x,y,errorv)

Name = ["$u$" "$v$"];
MaxError = [max(max(abs(erroru)));max(max(abs(errorv)))];
Norm  = [norm(erroru,"fro")/norm(zu,"fro");norm(errorv,"fro")/norm(zv,"fro")];
table(Norm , MaxError,'RowNames',Name)