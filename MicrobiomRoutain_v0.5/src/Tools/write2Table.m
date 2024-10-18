function write2Table(fileName,rowName,colName,Mat,Type)
if nargin<5
    Type = 'num';
end
f_id=fopen(strcat(fileName,'.txt'),'w');
M = length(rowName);
N = length(colName);
for i=1:N
    fprintf(f_id,'\t%s',colName{i});
end
fprintf(f_id,'\n');
switch Type
    case 'num'
        for i=1:M
            fprintf(f_id,'%s',rowName{i});
            for j=1:N
                fprintf(f_id,'\t%f',Mat(i,j));
            end
            fprintf(f_id,'\n');
        end
    case 'cell'
        for i=1:M
            fprintf(f_id,'%s',rowName{i});
            for j=1:N
                fprintf(f_id,'\t%s',Mat{i,j});
            end
            fprintf(f_id,'\n');
        end
end
end