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
if ~isfield(para,'thr')
    para.thr = 2;
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
    N_s = score(score<0);
    P_s = score(score>0);
    [~,N_id] = min(N_s);
    [~,P_id] = max(P_s);
    Enrich_lefse(score<0) = Neg_name(N_id);
     Enrich_lefse(score>0) = Pos_name(P_id);
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
    if ~isempty(score)
        plotSEL = abs(score)>=para.thr;
        plotLEfSe_pair(score(plotSEL),Enrich_lefse(plotSEL),tax_lefse(plotSEL), y_cmp)
        keyboard
        plotPDF(gcf,rpath);
    end
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
