%Tests + plots: stress

[x1, y1] = chebpts2(NX1, NY, [DOMX1(1) DOMX1(2) DOMY(1) DOMY(2)]);
[x2, y2] = chebpts2(NX2, NY, [DOMX2(1) DOMX2(2) DOMY(1) DOMY(2)]);

zSigmaX1=-(p./pi).*(-(2.*a.*y1.*(a.^2-x1.^2+y1.^2))./((y1.^2+(x1-a).^2).*(y1.^2+(x1+a).^2))+atan((a-x1)./y1)+atan((a+x1)./y1));
zSigmaY1=-(p./pi).*(+(2.*a.*y1.*(a.^2-x1.^2+y1.^2))./((y1.^2+(x1-a).^2).*(y1.^2+(x1+a).^2))+atan((a-x1)./y1)+atan((a+x1)./y1));
zShear1=-(4.*a.*p.*x1.*y1.^2)./(pi.*(y1.^2+(a+x1).^2).*(y1.^2+(a-x1).^2));

zSigmaX1(isnan(zSigmaX1))=-p;
zSigmaY1(isnan(zSigmaY1))=-p;
zShear1(isnan(zShear1))=0;

zSigmaX2=-(p./pi).*(-(2.*a.*y2.*(a.^2-x2.^2+y2.^2))./((y2.^2+(x2-a).^2).*(y2.^2+(x2+a).^2))+atan((a-x2)./y2)+atan((a+x2)./y2));
zSigmaY2=-(p./pi).*(+(2.*a.*y2.*(a.^2-x2.^2+y2.^2))./((y2.^2+(x2-a).^2).*(y2.^2+(x2+a).^2))+atan((a-x2)./y2)+atan((a+x2)./y2));
zShear2=-(4.*a.*p.*x2.*y2.^2)./(pi.*(y2.^2+(a+x2).^2).*(y2.^2+(a-x2).^2));

zSigmaX2(isnan(zSigmaX2))=0;
zSigmaY2(isnan(zSigmaY2))=0;
zShear2(isnan(zShear2))=0;

x=[x1,x2];
y=[y1,y2];
MATsigx=[MATsigx1,MATsigx2];
MATsigy=[MATsigy1,MATsigy2];
MATsh=[MATsh1,MATsh2];
zSigmaX=[zSigmaX1,zSigmaX2];
zSigmaY=[zSigmaY1,zSigmaY2];
zShear=[zShear1,zShear2];
%Error matrices
errorSigmaX = zSigmaX-MATsigx;
errorSigmaY = zSigmaY-MATsigy;
errorShear = zShear-MATsh;

subplot(2,3,1);
mesh(x,y,MATsigx)
subplot(2,3,2);
mesh(x,y,MATsigy)
subplot(2,3,3);
mesh(x,y,MATsh)
subplot(2,3,4);
mesh(x,y,errorSigmaX)
subplot(2,3,5);
mesh(x,y,errorSigmaY)
subplot(2,3,6);
mesh(x,y,errorShear)

Name = ["$\sigma_{x}$" "$\sigma_{y}$" "$\tau_{xy}$"];
MaxError = [max(max(abs(errorSigmaX)));max(max(abs(errorSigmaY)));max(max(abs(errorShear)))];
Norm  = [norm(errorSigmaX,"fro")/norm(zSigmaX,"fro");norm(errorSigmaY,"fro")/norm(zSigmaY,"fro");norm(errorShear,"fro")/norm(zShear,"fro")];
table(Norm , MaxError,'RowNames',Name)