function [FX,med,Q]=robust_scaling(X,p)
[d,N]=size(X);
if nargin<2
    p=0.05;
end
for i=1:d
    temp=X(i,:);
    if length(find(temp~=0))>0
Q(i,:) = quantile(X(i,temp~=0),[p 1-p],2);
    else
        Q(i,:) =[0 0];
    end
end
SQ=Q(:,2)-Q(:,1);
med=median(X,2);
FX=zeros(size(X));
for i=1:d
    st=std(X(i,(X(i,:)>=Q(i,1)&X(i,:)<=Q(i,2))));
    if st~=0
    FX(i,:)=(X(i,:)-med(i))/st;
    else
        FX(i,:)=(X(i,:)-med(i));
    end
end
FX(SQ==0,:)=0;
end