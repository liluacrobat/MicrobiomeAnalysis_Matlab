function A=ObserveOTU(D)
A=zeros(1,size(D,2));
for i=1:length(A)
    s=length(find(D(:,i)>0));
    A(i)=s;
end
end