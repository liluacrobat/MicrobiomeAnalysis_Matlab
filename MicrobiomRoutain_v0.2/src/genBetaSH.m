function genBetaSH(rpath_lefse,condaenv)
if nargin<2
    condaenv = '/Users/lli59/mothur';
end
fid = fopen(strcat(rpath_lefse,'.sh'),'w');
fprintf(fid,'#!/bin/sh\n')
fprintf(fid,'conda activate %s\n',condaenv)
t = strcat(rpath_lefse,'_Beta_div_tbl');
s = strsplit(t,'/');

fprintf(fid,'biom convert -i %s.txt -o %s.biom --to-json --table-type "OTU table"\n',s{end},s{end})
fprintf(fid,'mothur "#make.shared(biom=%s.biom)"\n',s{end})
fprintf(fid,'mothur "#dist.shared(shared=%s.shared, calc=thetayc-jclass-braycurtis, subsample=t)"\n',s{end})
end
