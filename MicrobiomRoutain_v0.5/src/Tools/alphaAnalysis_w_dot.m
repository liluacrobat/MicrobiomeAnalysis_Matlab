function alphaAnalysis_w_dot(fileName,measure_alpha,measure_alpha_title,groups_y,header)
plotBox(measure_alpha,groups_y,header);
ylabel(measure_alpha_title);
a=axis;
a(3) = max(a(3),min(measure_alpha(:))*0.8);

a(4)=a(4)*1.2;
axis(a);
keyboard
plotPDF(gcf,fileName);

fileName_anova = strcat(fileName,'_anova');
hg = header(groups_y);
[p_anova, t_anova,stats_anova] = anova1(measure_alpha,hg);
[c1,m,~,nms] = multcompare(stats_anova,'CriticalValueType','bonferroni'); 
[c2,~,~,~] = multcompare(stats_anova,'CriticalValueType','tukey-kramer'); 
alpha_pair_pval = zeros(max(groups_y));

c{1} = c1;
c{2} = c2;

for i=1:max(groups_y)
    for j=1:max(groups_y)
        [~,alpha_pair_pval(i,j)] = ttest2(measure_alpha(groups_y==i),measure_alpha(groups_y==j));
    end
end
writeMulCom(fileName,p_anova,t_anova,c,m,nms,hg,header,alpha_pair_pval);
end
function writeMulCom(filename,pval,t_anova,c,m,nms,hg,header,alpha_pair_pval)

y = zeros(1,length(nms));
for i=1:length(nms)
    for j=1:length(hg)
        if strcmp(nms{i},hg{j})==1
            y(i) = y(i)+1;
        end
    end
end
fid = fopen(strcat(filename,'.txt'),'w');
fprintf(fid,'One-way ANOVA F: %f\n',t_anova{2,5});
fprintf(fid,'One-way ANOVA p-value: %f\n',pval);

fprintf(fid,'\n');

fprintf(fid,'Groups\tGroup_ID\tMean\tStd\t# of samples\n');
for i=1:length(nms)
    fprintf(fid,'%d\t%s\t%f\t%f\t%d\n',i,nms{i},m(i,1),m(i,2),y(i));
end
fprintf(fid,'\n');
fprintf(fid,'Group-A\tGroup-B\tLower 95%% of diff (Group-A - Group-B)\tMean diff (Group-A - Group-B)\tHigher 95%% of diff (Group-A - Group-B)\n');
c1 = c{1};
c2 = c{2};

[m,~] = size(c1);
anova_p1 = ones(length(nms));
anova_p2 = ones(length(nms));

for i=1:m
    fprintf(fid,'%s\t%s\t%f\t%f\t%f\n',nms{c1(i,1)},nms{c1(i,2)},c1(i,3),c1(i,4),c1(i,5));
    anova_p1(c1(i,1),c1(i,2)) = c1(i,6);
    anova_p1(c1(i,2),c1(i,1)) = c1(i,6);

    anova_p2(c2(i,1),c2(i,2)) = c2(i,6);
    anova_p2(c2(i,2),c2(i,1)) = c2(i,6);
end
[hidx,~] = AlignID(header,nms);
anova_p1 = anova_p1(hidx,hidx);
anova_p2 = anova_p2(hidx,hidx);

nms_s = nms(hidx);
fprintf(fid,'\n');

fprintf(fid,'Pairwise ANOVA p-value (Bonferroni)\n');
for i=1:length(nms)
    fprintf(fid,'\t%s',nms_s{i});
end
fprintf(fid,'\n');
for i=1:length(nms)
    fprintf(fid,'%s',nms_s{i});
    for j=1:length(nms)
        fprintf(fid,'\t%f',anova_p1(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n');

fprintf(fid,'Pairwise ANOVA p-value (Tukey-Kramer)\n');
for i=1:length(nms)
    fprintf(fid,'\t%s',nms_s{i});
end
fprintf(fid,'\n');
for i=1:length(nms)
    fprintf(fid,'%s',nms_s{i});
    for j=1:length(nms)
        fprintf(fid,'\t%f',anova_p2(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n');

fprintf(fid,'Pairwise t-test p-value\n');
for i=1:length(nms)
    fprintf(fid,'\t%s',header{i});
end
fprintf(fid,'\n');
for i=1:length(nms)
    fprintf(fid,'%s',header{i});
    for j=1:length(nms)
        fprintf(fid,'\t%f',alpha_pair_pval(i,j));
    end
    fprintf(fid,'\n');
end
end
function plotBox(X,Y,Name,facecolor)
figure,
hold on;
wid = 0.3;
xidx = ~isnan(X);
X = X(xidx);
Y = Y(xidx);
if nargin<4
    facecolor = defaultColor(length(Name));
end
[Y,idx] = sort(Y);
X = X(idx);
boxplot(X,Y,'Colors',facecolor,'Widths',wid,'Symbol','w+');
h = findobj(gca,'Tag','Box');
for i=length(h):-1:1
    patch(get(h(i),'XData'),get(h(i),'YData'),facecolor(length(h)-i+1,:),'FaceAlpha',0.5,'EdgeColor',facecolor(length(h)-i+1,:));
end
% for i=length(h):-1:1
% %     patch(get(h(i),'XData'),get(h(i),'YData'),facecolor(length(h)-i+1,:),'FaceAlpha',0.5,'EdgeColor','none');
% end
h = boxplot(X,Y,'Colors',facecolor,'Widths',wid,'Symbol','w+');
UY = unique(Y);
for i=1:length(UY)
    tx = X(Y==UY(i));
    tt = (rand(length(tx),1)-0.5)*(wid*0.7)+i;
    plot(tt,tx,'o','MarkerFaceColor',facecolor(i,:),'MarkerSize',5,'MarkerEdgeColor',facecolor(i,:))
end
% set(h,{'linew'},{1.5});
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'linewidth',1);
set(gca,'FontSize',14);

if nargin>2
    set(gca,'xticklabels',Name);
end
if length(Name)<=3
    pbaspect([1 2 1])
else
    pbaspect([2 1 1])
end
a=axis;
a(1)=0.5;
a(2) = length(unique(Y))+0.5;
axis(a);
end
