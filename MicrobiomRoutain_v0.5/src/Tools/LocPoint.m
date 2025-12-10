function [R,Ri]=LocPoint(PCA,L,s)
[m,n]=size(L);
if nargin<3
    s=0.0001;
end
LMAX=L+s;
LMIN=L-s;
R=cell(1,n);
Ri=cell(1,n);
for i=1:n
    idx=ones(1,n);
    for j=1:m
        idx=idx&(PCA(j,:)>LMIN(j,i))&(PCA(j,:)<LMAX(j,i));
    end
    Ri{i}=find(idx==1);
    R{i}=PCA(:,idx==1);
end
end