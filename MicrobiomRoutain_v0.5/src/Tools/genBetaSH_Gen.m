function genBetaSH_Gen(path2tbl,condaenv)
if nargin<2
    condaenv = '/Users/lli59/mothur';
end
fid = fopen(('tmp_cal_beta.sh'),'w');
fprintf(fid,'#!/bin/sh\n')
fprintf(fid,'conda activate %s\n',condaenv)
s = strsplit(path2tbl,'/');

fprintf(fid,'biom convert -i %s.txt -o %s.biom --to-json --table-type "OTU table"\n',s{end},s{end})
fprintf(fid,'mothur "#make.shared(biom=%s.biom)"\n',s{end})
fprintf(fid,'mothur "#dist.shared(shared=%s.shared, calc=thetayc-jclass-braycurtis, subsample=t)"\n',s{end})
end
