function [FX,med,Q]=robust_scaling_norm(X,p)
[d,N]=size(X);
if nargin<2
    p=0.05;
end
for i=1:d
    if length(unique(X(i,:)))>1
        Q(i,:) = quantile(X(i,:),[p 1-p],2);
    else
        Q(i,:) =[0 0];
    end
end
SQ=Q(:,2)-Q(:,1);
med=median(X,2);
FX=zeros(size(X));
for i=1:d
    if SQ(i)~=0
        FX(i,:)=(X(i,:)-Q(i,1))/SQ(i);
    else
        FX(i,:)=0;
    end
end
FX(FX>1)=1;
FX(FX<0)=0;
end