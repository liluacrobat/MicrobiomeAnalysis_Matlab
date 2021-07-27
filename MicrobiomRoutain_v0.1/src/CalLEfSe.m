function [order,score,Enrich_lefse,tax_lefse,rpath] = CalLEfSe(x,y,para)
tax = para.tax;
y_legend = para.y_legend;
if ~iscell(y_legend)
    tmp = cell(size(y_legend));
    for i=1:length(tmp)
        tmp{i} = num2str(y_legend(i));
    end
    y_legend = tmp;
end
rpath = para.rpath;
rpath = strcat(rpath,'_CMP');
y_cmp = y_legend(unique(y));

u = unique(y);
for i=1:length(u)
    m_g(:,i) = mean(x(:,y==u(i)),2);
end
[~,Enrich_m] = max(m_g,[],2);
Enrich = y_legend(u(Enrich_m));

save('data3lefse_temp','x','y');
system('R CMD BATCH EsLEfSe.R');
load('lefse_ana_result_tmp','res_group');

order = mkList(res_group.Names);
score = res_group.scores;
Enrich_lefse = [];
tax_lefse = [];
fid = fopen(strcat(rpath,'.txt'),'w');
fprintf(fid,'Seq\tTax\tEnriched Group\tScore\n');
if ~isempty(res_group.Names)
    Enrich_lefse = Enrich(order);
    tax_lefse = tax(order);
    
    Neg_name = Enrich_lefse(score<0);
    Pos_name = Enrich_lefse(score>0);
    if length(Neg_name)>length(Pos_name)
        name_det = Neg_name;
        name_flag = 0;
    else
        name_det = Pos_name;
        name_flag = 1;
    end
    num = 0;
    for i=1:length(name_det)
        if strcmp(name_det{i},y_cmp{1})==1
            num=num+1;
        end
    end
    if num>length(name_det)
        if name_flag == 1
            name_as = y_cmp{1};
            name_as_n = y_cmp{2};
        else
            name_as = y_cmp{2};
            name_as_n = y_cmp{1};
        end
    else
        if name_flag == 1
            name_as = y_cmp{2};
            name_as_n = y_cmp{1};
        else
            name_as = y_cmp{1};
            name_as_n = y_cmp{2};
        end
    end
    for i=1:length(score)
        if score(i)>0
            Enrich_lefse{i} = name_as;
        else
            Enrich_lefse{i} = name_as_n;
        end
    end
    
    for i=1:length(Enrich_lefse)
        fprintf(fid,'%d\t%s\t%s\t%f\n',order(i),tax_lefse{i},Enrich_lefse{i},score(i));
    end
    
end
save(rpath,'order','score','Enrich_lefse','tax_lefse');
delete('data3lefse_temp.mat');
delete('lefse_ana_result_tmp.mat');
if isfield(para,'plot')==0
    para.plot = 1;
end
if para.plot~=0
    plotLEfSe(score,Enrich_lefse,tax_lefse, y_cmp)
    keyboard
    plotPDF(gcf,rpath);
end
end
function y = mkList(x)
if iscell(x)
    y = size(x(:));
    for i=1:length(x)
        y(i) = str2num(x{i});
    end
else
    if ischar(x)
        y = str2num(x);
    else
        error('Cannot identify the order');
    end
end
end
