function plotBarPoint(x,y,yl,FaceColor)
sel = ~isnan(x);
x=x(sel);
y=y(sel);
ori = 'vertical';
figure;
hold on
u = unique(y);
w=0.4;
for i=1:length(u)
    temp = x(y==u(i));
    s(i) = std(temp);
    m(i) = median(temp);
    s1(i) = quantile(temp,0.25);
    s2(i) = quantile(temp,0.75);
    
    %     errorbar(i,m(i),-s(i),s(i),ori,'Marker','none','Color','k','linewidth',1.5);
    plot([i-w*0.3,i+w*0.3],[s1(i) s1(i)],'-','Color','k','linewidth',1.5);
    plot([i-w*0.3,i+w*0.3],[s2(i) s2(i)],'-','Color','k','linewidth',1.5);
    plot([i-w*0.7,i+w*0.7],[m(i) m(i)],'-','Color','k','linewidth',1.5);
    plot([i,i],[s1(i) s2(i)],'-','Color','k','linewidth',1.5);
    
end
if nargin<4
    FaceColor = defaultColor;
end

rng default
for i=1:length(u)
    temp = x(y==u(i));
    temp1 = (temp-min(temp))/(max(temp)-min(temp));
    temp1 = temp1(:)';
    temp1(isnan(temp1)) = 0;
    r = (rand(1,length(temp))-0.5).*abs(1-temp1+mean(temp1))*w;
    plot(r+i,temp(:)', 'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor',FaceColor(i,:));
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
set(gca,'FontSize',12,'fontweight','bold');

% fig=gca;
% fig.xticks(fontsize=14)
% fig.YAxis.FontSize = 16;
% fig.XAxis.FontSize = 13 ;
% set(gca, 'XScale', 'log')
% xlabel('Relative abundance of \it{P. gingivalis}','FontSize',18);
set(gca,'FontSize',14);
pbaspect([1 1 1 ])
end