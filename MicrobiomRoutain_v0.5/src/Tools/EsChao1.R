rm(list=ls())
require(R.matlab)
require(fossil)
D<-readMat("Chao1_temp.mat");
X<-D$X;
# estimate Chao1 index
rich <- apply(X,2,chao1);
writeMat("Chao1_index.mat",rich=rich);
