function [groups_y, header] = genGroups(label,legend_ls)

[n_s,n_key] = size(label);
groups = sum(label.*repmat(10.^(0:n_key-1),n_s,1),2);
group_l = unique(groups);
groups_y = groups;
for i=1:length(group_l)
    groups_y(groups==group_l(i)) = i;
end
header = cell(1,length(group_l));

for i=1:length(group_l)
    gsel = gpL(group_l(i));
    t = legend_ls{1};
    if iscell(t)
        header{i} = t{gsel(1)};
    else
        header{i} = num2str(t(gsel(1)));
    end
    for j=2:length(gsel)
        t = legend_ls{j};
        header{i} = strcat(header{i},'-',t{gsel(j)});
    end
end
end
function y = gpL(x)
n = fix(log10(x))+1;
y = zeros(1,n);
for i=0:n-2
    y(i+1) = fix(x/10^i)-fix(x/10^(i+1))*10;
end
y(n) = fix(x/10^(n-1));
end