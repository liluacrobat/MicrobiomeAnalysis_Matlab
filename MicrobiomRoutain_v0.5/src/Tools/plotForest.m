function plotForest(fname,num)
ori = 'vertical';
figure;
hold on
FaceColor = defaultColor;
w=0.4;
if nargin<3
    cla=2;
end

% num = -num;
[~,idx] = sort(num,'ascend');
num = num(idx);

for i=1:length(num)
    s = abs(num(i)>0)+1;
    barh(i,num(i), 'FaceColor', FaceColor(s,:),'linewidth',1,'EdgeColor',FaceColor(s,:),'BarWidth',w);
end
a=axis;
a(1) = a(1)-0.2;
 a(2) = a(2)+0.2;
a(3) =0;
a(4) = length(num)+1;
axis(a);
xticks(-1:0.2:1);
xticklabels([1:-0.2:0 0.2:0.2:1]);
yticks([]);
xlabel('Importance Score');
% set(gca,'FontSize',12);
% set(gca,'FontSize',14);
pbaspect([1 1 1 ])
 plotPDF(gcf,fname);

end