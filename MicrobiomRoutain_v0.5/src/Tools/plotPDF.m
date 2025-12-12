function plotPDF(fig,filename)
savefig(fig, strcat(filename, '.fig'))
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,strcat(filename, '.pdf'),'-dpdf', '-r0');
% print(fig,strcat(filename, '.eps'),'-depsc', '-r600')
end