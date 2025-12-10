function R=OutportText(name,R,style)
if nargin<3
    style='num';
end
id=fopen(name,'w');
fprintf(id,' ');
for i=1:length(R.colname)
    fprintf(id,'\t%s',R.colname{i});
end
fprintf(id,'\n');
num=size(R.Data,2);
for i=1:length(R.rowname)
    fprintf(id,'%s',R.rowname{i});
    if strcmpi(style,'txt')
        for j=1:num
            fprintf(id,'\t%s',R.Data{i,j});
        end
    else
        for j=1:num
            fprintf(id,'\t%s',num2str(R.Data(i,j)));
        end
    end
    fprintf(id,'\n');
end
end