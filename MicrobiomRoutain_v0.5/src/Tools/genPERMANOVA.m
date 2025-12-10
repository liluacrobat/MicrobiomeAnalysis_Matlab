function genPERMANOVA(dm,Label)
fid = fopen('dissimilarity.txt','w');
N = size(dm,2);
for i=1:N
    fprintf(fid,'\t%d',i);
end
fprintf(fid,'\n');
for i=1:N
    fprintf(fid,'%d',i);
    for j=1:N
        fprintf(fid,'\t%f',dm(i,j));
    end
    fprintf(fid,'\n');
end
fid2 = fopen('dissimilarity_mapping.txt','w');
fprintf(fid2,'#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tCategory\tDescription\n');
for i=1:N
    fprintf(fid2,'%d\tAAAAAAAAAA\tAAAAAAAA\t%d\tNA\n',i,Label(i));
end
end