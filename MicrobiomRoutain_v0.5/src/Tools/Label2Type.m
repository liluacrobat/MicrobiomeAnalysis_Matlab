function R=Label2Type(L,U)
%% change the numeric label into diagnosis annotation
R=cell(size(L));
for i=1:length(R)
    R{i}=U{L(i)};
end
end