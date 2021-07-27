function R1=CalCLR(X)
R1=X;
if min(X(:))==0
    X=X+1;
end
for i=1:size(X,2)
    temp=X(:,i);
    temp2=temp(temp~=0);
    %     temp(temp==0)=0.1;
    m = geomean(temp2);
    R1(:,i)=log(temp)-log(m);
end
end