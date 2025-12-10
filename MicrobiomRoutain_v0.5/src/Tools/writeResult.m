function writeResult(map_name,OTU_ID,M,Diff,P,Q)
fid=fopen(map_name,'w');
fprintf(fid,'OTU ID\tMean(%%)\tDiff\tP\tQ\n');
for i=1:length(OTU_ID)
    fprintf(fid,'%s\t%f\t%f\t%f\t%f\n',OTU_ID{i},M(i),Diff(i),P(i),Q(i));
end
end







