function alphaAnalysis(fileName,measure_alpha,measure_alpha_title,groups_y,header)
plotBox(measure_alpha,groups_y,header);
ylabel(measure_alpha_title);
plotPDF(gcf,fileName);

fileName_anova = strcat(fileName,'_anova');
hg = header(groups_y);
[p_anova, t_anova,stats_anova] = anova1(measure_alpha,hg);
[c,m,~,nms] = multcompare(stats_anova,'CType','bonferroni');

alpha_pair_pval = zeros(max(groups_y));
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
fprintf(fid,'Group-A\tGroup-B\tLower 95%% of diff (Group-A - Group-B)\tMean diff (Group-A - Group-B)\tHigher 95%% of diff (Group-A - Group-B)\tp-value\n');
[m,~] = size(c);
anova_p = ones(length(nms));
for i=1:m
    fprintf(fid,'%s\t%s\t%f\t%f\t%f\t%f\n',nms{c(i,1)},nms{c(i,2)},c(i,3),c(i,4),c(i,5),c(i,6));
    anova_p(c(i,1),c(i,2)) = c(i,6);
    anova_p(c(i,2),c(i,1)) = c(i,6);
end
[hidx,~] = AlignID(header,nms);
anova_p = anova_p(hidx,hidx);
nms_s = nms(hidx);
fprintf(fid,'\n');

fprintf(fid,'Pairwise ANOVA p-value\n');
for i=1:length(nms)
fprintf(fid,'\t%s',nms_s{i});
end
fprintf(fid,'\n');
for i=1:length(nms)
    fprintf(fid,'%s',nms_s{i});
    for j=1:length(nms)
        fprintf(fid,'\t%f',anova_p(i,j));
    end
    fprintf(fid,'\n');
end

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
boxplot(X,Y,'Colors','k','Widths',wid,'Symbol','k+');
h = findobj(gca,'Tag','Box');
for i=length(h):-1:1
    patch(get(h(i),'XData'),get(h(i),'YData'),facecolor(i,:),'FaceAlpha',0.8,'EdgeColor','none');
end
h = boxplot(X,Y,'Colors','k','Widths',wid,'Symbol','k+');
set(h,{'linew'},{1.5});
lines = findobj(gcf,'type','line','Tag','Median');
set(lines,'linewidth',3);
set(gca,'FontSize',14);

if nargin>2
    set(gca,'xticklabels',Name);
end
if length(Name)<=3
    pbaspect([1 1 1])
else
    pbaspect([2 1 1])
end
a=axis;
a(1)=0.5;
a(2) = length(unique(Y))+0.5;
axis(a);
end
