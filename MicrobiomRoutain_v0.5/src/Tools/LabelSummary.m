function LabelSummary(Y,L)
UU=unique(Y);
U=L(UU);
for i=1:length(U)
    display([U{i} ': ' num2str(length(find(Y==UU(i))))]);
end
end