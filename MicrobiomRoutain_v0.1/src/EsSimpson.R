rm(list=ls())
require(R.matlab)
require(vegan)
D<-readMat("simpson_temp.mat");
X<-D$X;
# estimate Chao1 index
rich <-diversity(X, index = "simpson", MARGIN = 2)
writeMat("simpson_index.mat",rich=rich);
