% test + plot: stress

% SigmaX
z1=-(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
%sigmaY
z2=-(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
%SHEAR
z3=-(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

errorSigmaX = z1-MATsigx;
errorSigmaY = z2-MATsigy;
errorShear = z3-MATsh;

Name = {'Sigmax';'Sigmay';'Shear'};
MaxError = [max(max(abs(errorSigmaX)));max(max(abs(errorSigmaY)));max(max(abs(errorShear)))];
NormInf  = [norm(errorSigmaX,'Inf');norm(errorSigmaY,'Inf');norm(errorShear,'Inf')];

table(NormInf , MaxError,'RowNames',Name)

figure(5);
set(gcf,'units','points','position',[0,0,1000,600])
subplot(1,3,1)
mesh(x,y,errorSigmaX)
t=title('Error $\sigma_x$','FontSize',16);
set(t, 'interpreter', 'latex');
xlabel('$x$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$y$','interpreter','latex', 'FontWeight','bold','FontSize',16)
subplot(1,3,2)
mesh(x,y,errorSigmaY)
t=title('Error $\sigma_y$','FontSize',16);
set(t, 'interpreter', 'latex');
xlabel('$x$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$y$','interpreter','latex', 'FontWeight','bold','FontSize',16)
subplot(1,3,3)
mesh(x,y,errorShear)
t=title('Error $\tau_{xy}$','FontSize',16);
set(t, 'interpreter', 'latex');
xlabel('$x$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$y$','interpreter','latex', 'FontWeight','bold','FontSize',16)
%saveas(gcf, 'fig', 'epsc')
