function [L,U] = Type2Label(X)
%% change the numeric label into diagnosis annotation
U = unique(X);
uni = length(U);
n = length(X);
L = zeros(n,1)*nan;
for i=1:n
    for k=1:uni
        if strcmp(X{i},U{k})==1
            L(i) = k;
            break;
        end
    end
end
end