function plotBox(X,Y,Name,facecolor)
figure,
hold on;
wid = 0.3;
xidx = ~isnan(X);
X = X(xidx);
Y = Y(xidx);
if nargin<4
    facecolor = defaultColor;
end
[Y,idx] = sort(Y);
X = X(idx);
boxplot(X,Y,'Colors','k','Widths',wid,'Symbol','k+');
h = findobj(gca,'Tag','Box');
for i=length(h):-1:1
    patch(get(h(i),'XData'),get(h(i),'YData'),facecolor(mod(length(h)-i+1-1,7)+1,:),'FaceAlpha',0.8,'EdgeColor','none');
end
h = boxplot(X,Y,'Colors','k','Widths',wid,'Symbol','k+');
set(h,{'linew'},{1.5});
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'linewidth',3);
set(gca,'FontSize',14);

if nargin>2
    set(gca,'xticklabels',Name);
end
pbaspect([1 1 1])
a=axis;
a(1)=0.5;
a(2) = length(unique(Y))+0.5;
axis(a);
end
function colormap = defaultColor(n)
if nargin<1
    n=7;
end
if n<=7
colormap = [0, 0.4470, 0.7410
    0.8500, 0.3250, 0.0980
    0.9290, 0.6940, 0.1250
    0.4940, 0.1840, 0.5560
    0.4660, 0.6740, 0.1880
    0.3010, 0.7450, 0.9330
    0.6350, 0.0780, 0.1840];
else
    colormap = distinguishable_colors(n);
end
end