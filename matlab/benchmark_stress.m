% To test and plot stress fields

[x, y] = chebpts2(NX, NY, [DOMX(1) DOMX(2) DOMY(1) DOMY(2)]);
% SigmaX
z1=-(p./pi).*(-(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
% SigmaY
z2=-(p./pi).*(+(2.*a.*y.*(a.^2-x.^2+y.^2))./((y.^2+(x-a).^2).*(y.^2+(x+a).^2))+atan((a-x)./y)+atan((a+x)./y));
% Shear
z3=-(4.*a.*p.*x.*y.^2)./(pi.*(y.^2+(a+x).^2).*(y.^2+(a-x).^2));

% Error matrices
errorSigmaX = z1-MATsigx;
errorSigmaY = z2-MATsigy;
errorShear = z3-MATsh;

Name = ["$\sigma_{x}$" "$\sigma_{y}$" "$\tau_{xy}$"];
MaxError = [max(max(abs(errorSigmaX)));max(max(abs(errorSigmaY)));max(max(abs(errorShear)))];
Norm  = [norm(errorSigmaX,"fro")/norm(z1,"fro");norm(errorSigmaY,"fro")/norm(z2,"fro");norm(errorShear,"fro")/norm(z3,"fro")];
table(Norm , MaxError,'RowNames',Name)

% Save latex table
where=join(['..\Tab\',ScriptName,'.txt'],"");
file = fopen(where,'w');
fprintf(file, '%s\n%s\n%s\n','\begin{tabular}{c c c}','\toprule','\multicolumn{3}{l}{Parameters:} \\');
fprintf(file, '%s%d%s%d%s\n','$N_{x} \times N_{y}$& \multicolumn{2}{c}{',NX,'$\times$',NY,'}\\');
if NX+NY<81 % Condition number
    CondNumM=cond(M);
    fprintf(file, '%s%.3e%s\n','Condition no.& \multicolumn{2}{c}{',CondNumM,'}\\');
end
fprintf(file, '%s%.3f%s\n','Time& \multicolumn{2}{c}{',time,' s}\\');
fprintf(file, '%s\n%s\n%s\n','\midrule','Stress & Relative error & Max. absolute error \\','\midrule');
for i=1:3
    fprintf(file,'%s & %.3e & %.3e %s \n',Name(i),Norm(i),MaxError(i),'\\');
end
fprintf(file, '%s\n%s','\bottomrule','\end{tabular}');
fclose(file);

% Directory
newSubFolder=join(['../Fig/',ScriptName],'');
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end

% Save figure of error matrices
Z={errorSigmaX;errorSigmaY;errorShear};
index=["errorSigmaX" "errorSigmaY" "errorShear"];
for i = 1:3
   figure('Name',index(i))
   set(gcf,'units','centimeters','position',[0,0,14,10])
   mesh(x,y,cell2mat(Z(i)),'edgecolor', 'k');
   set(gca,'FontSize',10);
   xlabel('$x$','interpreter','latex', 'FontWeight','bold','FontSize',12)
   ylabel('$y$','interpreter','latex', 'FontWeight','bold','FontSize',12)
   where=join(["..\Fig\",ScriptName,"\Fig_",index(i),".pdf"],"");
   set(gca,'LooseInset',get(gca,'TightInset'));
   exportgraphics(gcf,where,'Resolution',300)
end

% Field plotter
Z={MATsigx;MATsigy;MATsh};
index=["sigmax" "sigmay" "shear"];
for i = 1:3
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
   for k=1:size(YTL,1)
       str2=sprintf('%1.0e',str2double(cell2mat(YTL(k))));
       if str2(1)=='-'
           YTL(k)={str2(1:2)};
       else
           YTL(k)={str2(1)};
       end
   end
   set(col,'YTickLabel','');
   set(col,'YTickLabel',YTL);
   if abs(DOMX(1)-DOMX(2))<c
       text(c+5,-6,join(['\times10^{',exponent,'}'],''))
   else
       text(c+15,-6,join(['\times10^{',exponent,'}'],''))
   end
   % Contour plot
   if Norm(i)<0.55
       contour(x,y,cell2mat(Z(i)),20,'k','LineWidth',1.0)
   elseif (Norm(i)>=0.4) && (Norm(i)<0.8)
       contour(x,y,cell2mat(Z(i)),10,'k','LineWidth',1.0)
   else
       contour(x,y,cell2mat(Z(i)),5,'k','LineWidth',1.0)
   end
   where=join(["..\Fig\",ScriptName,"\Fig_",index(i),".pdf"],"");
   set(gca, 'XAxisLocation', 'top')
   set(gca, 'YDir','reverse')
   set(gca,'LooseInset',get(gca,'TightInset'));
   hold off
   exportgraphics(gcf,where,'ContentType','image','Resolution',300)
end