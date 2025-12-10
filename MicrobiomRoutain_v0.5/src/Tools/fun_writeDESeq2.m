function fun_writeDESeq2(filename,T,ID,idname,gene_length)
fid = fopen(strcat(filename,'.txt'),'w');
fprintf(fid,'ID\tGeneName\tlog2FC\tFDR\tmeanExpAve\tgene_length\n');
for i=1:size(T,1)
    fprintf(fid,'%s',ID{i});
    fprintf(fid,'\t%s',idname{i});
    for j=1:size(T,2)
    fprintf(fid,'\t%f',T(i,j));
    end
    fprintf(fid,'\t%f',gene_length(i));
    fprintf(fid,'\n');
end

end