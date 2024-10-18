function R=CalRel(X,p)
if nargin<3
    R=X./repmat(sum(X,1),size(X,1),1);
else
    Q=zeros(1,size(X,2));
    for i=1:length(Q)
        temp=X(:,i);
        S=cumsum(sort(temp,'ascend'));
        Q(i)=quantile(S,p);
        if  Q(i)==0
            Q(i)=eps;
        end
    end
    R=X./repmat(Q,size(X,1),1);
end
end