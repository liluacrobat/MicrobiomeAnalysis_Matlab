function x = rmStr(x,tar)
for i=1:length(x)
    x{i} = strrep(x{i},tar,' ');
end
end