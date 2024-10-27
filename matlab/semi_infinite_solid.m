% To plot the stress and displacement analytical solutions

clear all

% Spatial dimensions
% Rectangle [-c,c]x[0,b]
b = 100;
c = 100;
a = 10;
% Acting pressure
p = 1e+06;
% Material constants
rat = 0.49;
E = 10e+09;
% rat = 0.25;
% E = 200e+09; 
lambda = E*rat/((1+rat)*(1-2*rat));
mu = E/(2*(1+rat));

% Steps
step1 = c/500;
step2 = b/500;

% Domain
[x, y] = meshgrid(-c:step1:c,0:step2:b);

% Analytical solutions
% SigmaX
z1=-(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
% SigmaY
z2=-(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
% ShearXY
z3=-(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));
% Max-shear
z4 = sqrt((z1./2-z2./2).^2+z3.^2);
% u
z5 = (1+rat)./(pi.*E)*p.*(-(2.*rat-1).*(a-x).*atan((a-x)./y)+(2.*rat-1).*(a+x).*atan((a+x)./y)+(rat-1).*y.*log((y.^2+(a-x).^2)./(y.^2+(a+x).^2)));
% Zero level
zero_level=@(x,y)  (2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2));
%v
% Options: zero_level(c,b) or zero_level(0,b)
z6 = (1+rat)./(pi.*E)*p.*((2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2))-zero_level(c,b));


% Plotter
Z={z1;z2;z3;z4;z5;z6};
index=["sigma_x" "sigma_y" "shear" "maxshear" "u" "v"];
for i = 1:6
   figure('Name',index(i))
   set(gcf,'units','centimeters','position',[0,0,14,10])
   s = pcolor(x,y,cell2mat(Z(i)));
   hold on
   set(s,'FaceColor', 'interp');
   set(s, 'EdgeColor', 'none');
   set(gca,'FontSize',10);
   xlabel('$x$','interpreter','latex', 'FontWeight','bold','FontSize',12)
   ylabel('$y$','interpreter','latex', 'FontWeight','bold','FontSize',12)
   col=colorbar;
   colormap(gca,parula)
   % Colorbar Exponent
   zmax=max(max(cell2mat(Z(i))));
   zmin=min(min(cell2mat(Z(i))));
   str=sprintf('%0.1e',(abs(zmax)+abs(zmin))/2);
   exponent=str(5:end);
   YTL = get(col,'YTickLabel');
   set(col,'YTickLabel','');
   set(col,'YTickLabel',YTL);
   text(c+15,-6,join(['\times10^{',exponent,'}'],''))
   % Contour plot
   contour(x,y,cell2mat(Z(i)),30,'k','LineWidth',0.7)
   where=join(["..\Fig\Fig_uniform_",index(i),".pdf"],"");
   set(gca, 'XAxisLocation', 'top')
   set(gca, 'YDir','reverse')
   set(gca,'LooseInset',get(gca,'TightInset'));
   hold off
   exportgraphics(gcf,where,'ContentType','image','Resolution',300)
end