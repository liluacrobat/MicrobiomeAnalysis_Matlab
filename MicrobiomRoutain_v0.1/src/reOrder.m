function [header,t] = reOrder(idx,header,y)
header = header(idx);
t = y;
for i=1:length(idx)
    t(y==idx(i)) = i;
end
end