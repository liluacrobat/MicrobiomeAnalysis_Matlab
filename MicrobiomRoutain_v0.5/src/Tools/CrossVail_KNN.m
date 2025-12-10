function [Fmacro,Fmicro,ACC]=CrossVail_KNN(X,label_or,knum)
[label,FeatureMapping]=Num2Label(label_or); 
fold_id=selectFOLD(label);
Test_Result=zeros(1,10);
w_history=zeros(size(X,1),10);
Para4Logo.plotfigure=0;
nk=length(unique(label));
for k=1:10
    display(['Forld ' num2str(k)]);
    train_patterns=X(:,fold_id~=k);
    train_targets=label(fold_id~=k);
    test_patterns=X(:,fold_id==k);
    test_targets=label(fold_id==k);
    Prediction= knnclassify(test_patterns', train_patterns', train_targets, knum);
    Test_Error = length(find(Prediction(:)~=test_targets(:)))/length(test_targets);
    Test_Result(k)=Test_Error;
    [Fmacro(k),Fmicro(k)]=F_Measure(Prediction,test_targets,nk);
end
R=(Test_Result*100);
ACC=100-R;
ACC=mean(ACC);
Fmacro=mean(Fmacro);
Fmicro=mean(Fmicro);
end
function l=selectLABEL(idx,s)
% select the labels included in set s
l=zeros(size(idx));
for i=1:length(s)
    l(idx==s(i))=1;
end
end
function Com=selectFOLD(label)
% random select samples for 10 fold test
idx=1:length(label);
rng default
nQ=length(unique(label));
l_G=zeros(1,nQ);
perm_G=cell(1,nQ);
step_G=zeros(1,nQ);
idx_G=cell(1,nQ);
for i=1:nQ
    l_G(i)=length(find(label==i));
    perm_G{i}=randperm(l_G(i));
    step_G(i)=round(l_G(i)/10);
    idx_G{i}=idx(label==i);
end
Com=ones(length(label),1)*10;
for i=1:9
    for j=1:length(unique(label))
        idx_B=idx_G{j};
        perm_B=perm_G{j};
        l_B=l_G(j);
        Com(idx_B(perm_B(round(l_B/10*(i-1))+1:round(l_B/10*i))))=i;
    end
end
end
function [Fmacro,Fmicro]=F_Measure(P,G,K)
TP=zeros(1,K);
FP=zeros(1,K);
FN=zeros(1,K);
N=length(G);
for i=1:K
    Tnum=length(find(G==i));
    Nnum=N-Tnum;
    temp=G(P==i);
    TPnum=length(find(temp==i));
    FPnum=length(find(temp~=i));
    temp=G(P~=i);
    FNnum=length(find(temp==i));
    TP(i)=TPnum/Tnum;
    FP(i)=FPnum/Nnum;
    FN(i)=FNnum/Tnum;
end
Pmicro=sum(TP)/sum(TP+FP);
Rmicro=sum(TP)/sum(TP+FN);
sum_PR_mic=Pmicro+Rmicro;
if sum_PR_mic~=0
    Fmicro=2*Pmicro*Rmicro/(Pmicro+Rmicro);
else
    Fmicro=2*Pmicro*Rmicro/(sum_PR_mic+eps);
end
Pmacro=TP./(TP+FP);
Pmacro(isnan(Pmacro))=0;
Pmacro=sum(Pmacro)/K;
Rmacro=TP./(TP+FN);
Rmacro(isnan(Rmacro))=0;
Rmacro=sum(Rmacro)/K;
sum_PR_mac=Pmacro+Rmacro;
if sum_PR_mac~=0
    Fmacro=2*Pmacro*Rmacro/(Pmacro+Rmacro);
else
    Fmacro=2*Pmacro*Rmacro/(sum_PR_mac+eps);
end
end
function [Label,U]=Num2Label(L)
Label=zeros(size(L,1),1);
U=unique(L);
for i=1:length(U)
    Label(L==U(i))=i;
end
end