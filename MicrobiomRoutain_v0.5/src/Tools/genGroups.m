function [groups_y, header] = genGroups(label,legend_ls)

[n_s,n_key] = size(label);

for i=1:n_key
    if ~iscell(legend_ls{i})
        tmp = legend_ls{i};
        tmpC = cell(size(tmp));
        for j=1:length(tmp)
            tmpC{j} =  num2str(tmp(j));
        end
    end
end
new_label = cell(n_s,1);
for i=1:n_s
    T1 = legend_ls{1};
    new_label{i}=T1{label(i,1)};
    for j=2:n_key
         T1 = legend_ls{j};
    new_label{i}=strcat(new_label{i},'-',T1{label(i,j)});
    end
end
[header,~,groups_y] = unique(new_label);
end