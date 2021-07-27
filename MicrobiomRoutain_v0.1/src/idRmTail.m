function y  = idRmTail(x)
y =cell(size(x));
for i=1:length(x)
    temp = x{i};
    s = strsplit(temp,'_S');
    s1 = temp(1:end-length(s{end})-2);
    y{i} = s1;
end
end