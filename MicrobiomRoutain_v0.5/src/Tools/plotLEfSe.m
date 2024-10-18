function plotLEfSe(score,Enrich_lefse,tax_lefse, legend_ls,FaceColor)
if nargin<5
    FaceColor = defaultColor(length(legend_ls));
end
[idx1, ~] = AlignID(Enrich_lefse, legend_ls);
w=0.4;
bx = 1:length(idx1);
bx = flip(bx);
u = unique(idx1);
figure, hold on;
for i=1:length(u)
    barh(bx(idx1==u(i)),score(idx1==u(i)),'FaceColor', FaceColor(i,:),'linewidth',1,'EdgeColor',FaceColor(i,:),'BarWidth',w);
end
yticks(1:length(idx1));
yticklabels(flip(tax_lefse));
xlabel('LDA SCORE (log 10)');
% ax = axis;
ax(1) = min(score)*1.2;
if length(unique(Enrich_lefse))==1
    ax(1) = 0;
end
ax(2) = max(score)*1.2;
ax(3) = 0;
ax(4) = length(idx1)+1;
if ax(1)==ax(2)
    ax(1) = 0;
end
axis(ax);
xt = xticks;
xticklabels(xt);
box on
legend(legend_ls(u));
set(gca,'FontSize',18);

if length(xt)<4
    xl = fix(max(abs(score)));
    xticks(-xl:xl);
    xticklabels((-xl:xl));
end
if length(score)>3
    if length(score)>30
        pbaspect([1 2 1]);
    else
        pbaspect([1 1.5 1]);
    end
else
    pbaspect([3 1 1]);
end


end
function [idx12,idx21]=AlignID(ID1,ID2)
% 2->1
n1=length(ID1);
n2=length(ID2);
idx12=zeros(n1,1);
idx21=zeros(n2,1);
for i=1:n1
    s1=strtrim(ID1{i});
    for j=1:n2
        s2=strtrim(ID2{j});
        if strcmpi(s1,s2)
            idx12(i)=j;
            break;
        end
    end
end
for i=1:n2
    s2=strtrim(ID2{i});
    for j=1:n1
        s1=strtrim(ID1{j});
        if strcmpi(s1,s2)
            idx21(i)=j;
            break;
        end
    end
end
end