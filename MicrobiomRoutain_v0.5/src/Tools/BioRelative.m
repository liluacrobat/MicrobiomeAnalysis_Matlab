function X=BioRelative(X)
X=X./repmat(sum(X,1),size(X,1),1);
end