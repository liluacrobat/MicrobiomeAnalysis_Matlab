function x = catLabel(meta,idx)
% Categorize string labels into 1,2,3,...
x.raw = table2array(meta(:,idx));
[x.legend, ~, x.y] = unique(x.raw);
end