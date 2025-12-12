function R1=CalCLR(X)
R1=X;
X = CalRel(X);
if min(X(:))==0
    tmp = X(:);
    X=min(tmp(tmp>0))/2+X;
end
for i=1:size(X,2)
    temp=X(:,i);
    temp2=temp(temp~=0);
    m = geomean(temp2);
    R1(:,i)=log2(temp)-log2(m);
end
end