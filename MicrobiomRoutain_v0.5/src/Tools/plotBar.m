function plotBar(x,y,yl,err,FaceColor)
sel = ~isnan(x);
x=x(:);
y=y(:);
x=x(sel);
y=y(sel);
ori = 'vertical';
figure;
hold on
u = unique(y);
if nargin<4
    err=0;
end
for i=1:length(u)
    temp = x(y==u(i));
    s(i) = std(temp);
    m(i) = mean(temp);
    if err==1
        errorbar(i,m(i),0,s(i),ori,'Marker','none','Color','k','linewidth',1);
    end
end
if nargin<5
    FaceColor = defaultColor;
end
w=0.4;
for i=1:length(u)
    bar(i, m(i), 'FaceColor', FaceColor(i,:),'BarWidth',w);
end
xticks(1:length(u));
xticklabels(yl(u));
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
pbaspect([1 1 1 ])
ax = axis;
ax(1) = 0;
ax(2) = length(u)+1;
axis(ax)
end