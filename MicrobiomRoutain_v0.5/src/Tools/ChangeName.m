function [New, idx]= ChangeName(OTU1,OTU2)
R1 = ExtractOTUid(OTU1);
R2 = ExtractOTUid(OTU2);
[idx,~] = AlignID(R1,R2);
if isempty(find(idx==0, 1))
    New = OTU2(idx);
else
    disp('Mismatch!!');
end
end
function R = ExtractOTUid(OTU)
for i=1:length(OTU)
    temp = OTU{i};
    s =strsplit(temp,'OTU_');
    R{i} = s{2};
end
end