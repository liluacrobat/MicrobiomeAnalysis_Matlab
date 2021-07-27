function [order,score,Enrich_lefse,tax_lefse,rpath] = genMultiLEfSe(x,y,para)
tax = para.tax;
y_legend = para.y_legend;
rpath = para.rpath;
rpath = strcat(rpath,'_CMP');
y_cmp = y_legend(unique(y));
% for i=1:length(y_cmp)
%     if i==1
%         rpath = strcat(rpath,y_cmp{i});
%     else
%         rpath = strcat(rpath,'-',y_cmp{i});
%     end
% end

filename = strcat(rpath,'_4Galaxy.txt');
fid = fopen(filename,'w');
fprintf(fid,'Sample_ID');
[m,n] = size(x);

UY = unique(y);

x4g = [];
y4g = [];
y = y(:);
for i=1:length(UY)
    x4g = [x4g x(:,y==UY(i))];
    y4g = [y4g;y(y==UY(i))];
end

for i=1:n
    fprintf(fid,'\tS%d',i);
end
fprintf(fid,'\n');

fprintf(fid,'Classes');
for i=1:n
    fprintf(fid,'\t%d',y4g(i));
end
fprintf(fid,'\n');

for i=1:m
    fprintf(fid,'OTU%d',i);
    for j=1:n
        fprintf(fid,'\t%d',x4g(i,j));
    end
    fprintf(fid,'\n');
end
keyboard
Galaxy_tbl = readtable(strcat(rpath,'_4Galaxy','.lefse_internal_res'),'delimiter','\t','FileType','delimitedtext');

OTU = table2array(Galaxy_tbl(:,1));
group_or = table2array(Galaxy_tbl(:,3));
score_or = table2array(Galaxy_tbl(:,4));
pval_or = table2array(Galaxy_tbl(:,5));
for i=1:length(OTU)
    order_or(i,1) = str2num(strrep(OTU{i},'OTU',''));
end
sel = ~isnan(score_or);

score_t = score_or(sel);
group_t = group_or(sel);
order_t = order_or(sel);

score =[];
order =[];
group = [];
UG = unique(group_t);
for i=1:length(UG)
    st = score_t(group_t==UG(i));
    gt = group_t(group_t==UG(i));
    ot = order_t(group_t==UG(i));
    
    [~,sort_idx] = sort(st,'descend');
    
    score =[score; st(sort_idx)];
    order =[order; ot(sort_idx)];
    group = [group; gt(sort_idx)];
    
end

Enrich_lefse = para.y_legend(group);
tax_lefse = tax(order);
save(rpath,'order','score','Enrich_lefse','tax_lefse');

plotLEfSe(score,Enrich_lefse,tax_lefse, para.y_legend)
plotPDF(gcf,rpath);
end
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
pbaspect([1 1.5 1]);
% ax = axis;
ax(1) = 0;
ax(2) = max(score)*1.2;
ax(3) = 0;
ax(4) = length(idx1)+1;
axis(ax);
xt = xticks;
if length(xt)<4
    xl = fix(max(abs(score)));
    xticks(-xl:xl);
end
box on
legend(legend_ls(u),'location','southeast');
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