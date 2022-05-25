function [y,h] = changeHeaderOrder(group,header,order)
h = header(order);
y=zeros(size(group))*nan;
for i=1:length(order)
    y(group==order(i)) = i;
end
if sum(isnan(y))>0 || length(unique(order))~=length(order)
    h = header;
    y = group;
    disp('Incorrect input');
end

end