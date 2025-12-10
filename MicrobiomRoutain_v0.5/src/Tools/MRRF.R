rm(list=ls())
require(R.matlab)
library(RRF);

D<-readMat("RRF_temp.mat");

set.seed(D$seed)

X<-D$X;
class <- D$Y
gamma <- D$gama
ntree <- D$ntree
#ordinary random forest.
rf <- RRF(X,as.factor(class), flagReg = 0)
impRF <- rf$importance
impRF <- impRF[,"MeanDecreaseGini"]

#guided regularized random forest
imp <- impRF/(max(impRF))#normalize the importance score
coefReg <- (1-gamma)+gamma*imp #weighted average
grrf <- RRF(X,as.factor(class),coefReg=coefReg, importance=TRUE, nPerm=1, mtry= D$p,flagReg=1,ntree=ntree)
feaSet <-grrf$feaSet
meanDeAcc<- grrf$importance[,'MeanDecreaseAccuracy']
meanDeAccSD <- grrf$importanceSD[,'MeanDecreaseAccuracy']
DeGini <- grrf$importance[,'MeanDecreaseGini']

if (D$preflag==1){
    TestX = D$TestSet
    response <- predict(grrf, TestX)
    response <- c(response)
    score <- predict(grrf, TestX,type="vote",norm.votes=FALSE)
    writeMat("RRF_result.mat",feaSet=feaSet,meanDeAcc=meanDeAcc,meanDeAccSD=meanDeAccSD,
    DeGini=DeGini,response=response,score=score);
}else{
    writeMat("RRF_result.mat",feaSet=feaSet,meanDeAcc=meanDeAcc,meanDeAccSD=meanDeAccSD,DeGini=DeGini);
}

save(grrf, file = D$filename)
