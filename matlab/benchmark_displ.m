% To test and plot displacement fields

[x, y] = chebpts2(NX, NY, [DOMX(1) DOMX(2) DOMY(1) DOMY(2)]);
% u
z5 = (1+rat)./(pi.*E)*p.*(-(2.*rat-1).*(a-x).*atan((a-x)./y)+(2.*rat-1).*(a+x).*atan((a+x)./y)+(rat-1).*y.*log((y.^2+(a-x).^2)./(y.^2+(a+x).^2)));
% Zero level
zero_level=@(x,y)  (2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2));
% v
z6 = (1+rat)./(pi.*E)*p.*((2.*rat-1).*y.*(atan((a-x)./y)+atan((a+x)./y))+(rat-1).*((a-x).*log(y.^2+(a-x).^2)+(a+x).*log(y.^2+(a+x).^2))-zero_level(c,b));

% Error matrices
erru=z5-MATu;
errv=z6-MATv;

Name = ["$u$" "$v$"];
Norm = [norm(erru,"fro")/norm(z5,"fro");norm(errv,"fro")/norm(z6,"fro")];
MaxError = [max(max(abs(erru)));max(max(abs(errv)))];
table(Norm,MaxError, 'RowNames',Name)

% Save a latex table
where=join(['..\Tab\',ScriptName,'.txt'],"");
file = fopen(where,'w');
fprintf(file, '%s\n%s\n%s\n','\begin{tabular}{c c c}','\toprule','\multicolumn{3}{l}{Parameters:} \\');
fprintf(file, '%s%d%s%d%s\n','$N_{x} \times N_{y}$& \multicolumn{2}{c}{',NX,'$\times$',NY,'}\\');
if NX+NY<161 % Condition number
    CondNumM=cond(M);
    fprintf(file, '%s%.3e%s\n','Condition no.& \multicolumn{2}{c}{',CondNumM,'}\\');
end
fprintf(file, '%s%.3f%s\n','Time& \multicolumn{2}{c}{',time,' s}\\');
fprintf(file, '%s\n%s\n%s\n','\midrule','Displacement & Relative error & Max. absolute error \\','\midrule');
for i=1:2 
    fprintf(file,'%s & %.3e & %.3e %s \n',Name(i),Norm(i),MaxError(i),'\\');
end
fprintf(file, '%s\n%s','\bottomrule','\end{tabular}');
fclose(file);

% Directory
newSubFolder=join(['../Fig/',ScriptName],'');
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end

% Save figures of error matrices
Z={erru;errv};
index=["erru" "errv"];
for i = 1:2
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
Z={MATu;MATv};
index=["u" "v"];
for i = 1:2
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
   text(c+5,-6,join(['\times10^{',exponent,'}'],''))
   % Contour plot
   contour(x,y,cell2mat(Z(i)),20,'k','LineWidth',1.0)
   where=join(["..\Fig\",ScriptName,"\Fig_",index(i),".pdf"],"");
   set(gca, 'XAxisLocation', 'top')
   set(gca, 'YDir','reverse')
   set(gca,'LooseInset',get(gca,'TightInset'));
   hold off
   exportgraphics(gcf,where,'ContentType','image','Resolution',300)
end
