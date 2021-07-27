function [head,tax] = rmHead(x)
s = strsplit(x,'__');
head = s{1};
tax = s{2};
end