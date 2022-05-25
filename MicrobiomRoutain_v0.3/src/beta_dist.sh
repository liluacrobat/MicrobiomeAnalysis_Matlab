#!/bin/sh
module use /projects/academic/pidiazmo/projectmodules
module load mothur/1.44.3
module load qiime2

biom convert -i DM_GG_table_ready.txt -o DM_GG.16S.L7.biom --to-json --table-type "OTU table"
mothur "#make.shared(biom=DM_GG.16S.L7.biom)"
mothur "#dist.shared(shared=DM_GG.16S.L7.shared, calc=thetayc-jclass-braycurtis, subsample=t)"

