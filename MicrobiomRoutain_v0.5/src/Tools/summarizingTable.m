function RS = summarizingTable(savename,Sample_ID,Data,tax_b,level,ppFlag)
if nargin<5
    level=0;
    ppFlag = 0;
end
RS = cell(1,length(level));
for i=1:length(level)
    sl = level(i);
    [mtx,tax_a] = fun_mergeLevel(Data,tax_b,sl);
    R.sample = Sample_ID;
    R.mtx = mtx;
    R.tax = tax_a;
    R.level = sl;
    if ppFlag==1
        if level==0
            writeOTUtable(strcat(savename,'_level_strain'),R.tax,R.sample,R.mtx);
        else
            writeOTUtable(strcat(savename,'_level_',num2str(sl)),R.tax,R.sample,R.mtx);
        end
    end
    RS{i} = R;
end
end
function [mtx,tax_a] = fun_mergeLevel(Data,tax_b,level)
% Data: data matrix
% tax_b: taxonomy before merging

if nargin<3
    level=0;
end
if level==0
    tax_a=unique(tax_b);
    mtx = fun_merging(Data,tax_b,tax_a);
else
    tax_a = tax_b;
    for i=1:length(tax_b)
        s = strsplit(tax_b{i},';');
        tax_a{i} = strjoin(s(1:level),';');
    end
    tax_a = unique(tax_a);
    mtx = fun_merging(Data,tax_b,tax_a);
end
end
function mtx = fun_merging(dd,tax_b,tax_a)
mtx = zeros(length(tax_a),size(dd,2));
parfor i=1:length(tax_a)
    ia = ismember(tax_b,tax_a{i})
    mtx(i,:) = sum(dd(ia,:),1);
end
end
function writeOTUtable(map_name,OTU_ID,Sample_ID,Data,tax)
Data = full(Data);
fid=fopen([map_name '.txt'],'w');
fprintf(fid,'#OTU_ID');
for i=1:length(Sample_ID)
    fprintf(fid,'\t%s',Sample_ID{i});
end
if nargin==5
    fprintf(fid,'\ttaxonomy');
end
fprintf(fid,'\n');
for i=1:length(OTU_ID)
    fprintf(fid,'%s',OTU_ID{i});
    for j=1:length(Sample_ID)
        fprintf(fid,'\t%f',Data(i,j));
    end
    if nargin==5
        fprintf(fid,'\t%s',tax{i});
    end
    fprintf(fid,'\n');
end
end

