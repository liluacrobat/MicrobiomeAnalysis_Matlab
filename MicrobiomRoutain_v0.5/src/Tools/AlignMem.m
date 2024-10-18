function [idx12,idx21]=AlignMem(ID1,ID2)
% 2->1
n1=length(ID1);
n2=length(ID2);
idx12=zeros(n1,1);
idx21=zeros(n2,1);
[la1,lc1] = ismember(ID1,ID2);
idx12(la1) = lc1(la1);
[la2,lc2] = ismember(ID2,ID1);
idx21(la2) = lc2(la2);
end