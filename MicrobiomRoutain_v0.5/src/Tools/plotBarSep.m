function plotBarSep(x,y,s,FaceColor)
sel = ~isnan(x);
x=x(sel);
y=y(sel);
ori = 'vertical';
figure;
hold on
u = unique(y);

if nargin<5
    FaceColor = defaultColor;
end
w=0.5;
n=0;
for i=1:length(u)
    tx = x(y==u(i));
    ty = n+(1:length(tx));
    n = ty(end);
    name(ty) = s(y==u(i));
    bar(ty, tx, 'FaceColor', FaceColor(i,:),'linewidth',1,'BarWidth',w);
end
xticks(1:n);
xticklabels(name);
% a=axis;
% a(3) = 0.3;
% a(4) = 2.7;
% a(1) = 10^-7;
% a(2) = 10^-1;
% axis(a);

box on
set(gca,'FontSize',12);

% fig=gca;
% fig.xticks(fontsize=14)
% fig.YAxis.FontSize = 16;
% fig.XAxis.FontSize = 13 ;
% set(gca, 'XScale', 'log')
% xlabel('Relative abundance of \it{P. gingivalis}','FontSize',18);
% set(gca,'FontSize',14);
pbaspect([3 1 1 ])
end