function [head,tax] = rmHead(x)
s = strsplit(x,'__');
if length(s)>1
    head = s{1};
    tax = s{2};
else
    head = '';
    tax = s{1};
end
end