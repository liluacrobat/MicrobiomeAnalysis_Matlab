function [pp,pval,rowname] = fun_ConPERMANOVA_rnd(mm,y,vadj,method,subj)
measure = mm;
group = [y vadj];
nvar = size(group,2);
save tmp4permanova measure group subj
fid = fopen('tmp_confoundPERMANOVA.R','w');
fprintf(fid,'rm(list=ls())\n')
fprintf(fid,'require(R.matlab)\n')
fprintf(fid,'library(vegan)\n')
fprintf(fid,'set.seed(123) \n')
fprintf(fid,'D<-readMat("tmp4permanova.mat")\n')
fprintf(fid,'Y <- data.frame(D$measure)\n')
fprintf(fid,'samp <- data.frame(D$group)\n')
fprintf(fid,'ps <- adonis2(Y ~ ')
fprintf(fid,'samp$X1');
for i=2:nvar
    fprintf(fid,'+samp$X%d',i);
end
fprintf(fid,', method = "%s",strata=D$subj,permutations = 1000)\n',method);
fprintf(fid,'pval <- ps$`Pr(>F)`\n');
fprintf(fid,'rn <- rownames(ps)\n');
fprintf(fid,'writeMat("confoundingPERMANOVA_result_tmp.mat",pval=pval,rowname=rn)')

system('R CMD BATCH tmp_confoundPERMANOVA.R');
load('confoundingPERMANOVA_result_tmp.mat','pval','rowname');
pp = pval(1);
delete('tmp4permanova.mat');
delete('confoundingPERMANOVA_result_tmp.mat');
end
