function R=transform_rankdata(X)
R=X;
for i=1:size(X,1)
    R(i,:)=Rankings(X(i,:));
end
end