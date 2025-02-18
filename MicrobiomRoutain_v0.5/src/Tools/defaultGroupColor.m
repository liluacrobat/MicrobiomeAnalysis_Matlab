function colorGroups = defaultGroupColor(n)
k = length(n);
n(k+1:7) = 7;
colorGroups{1} =  defaultColor(n(1));
% colorGroups{1} =  cbrewer('qual', 'Set1',n(1));
colorGroups{2} =  cbrewer('qual', 'Set2',n(2));
colorGroups{3} =  cbrewer('qual', 'Set3',n(3));
colorGroups{4} =  cbrewer('qual', 'Paired',n(4));
colorGroups{5} =  cbrewer('qual', 'Accent',n(5));
colorGroups{6} =  cbrewer('qual', 'Dark2',n(6));
colorGroups{7} = [0, 0.4470, 0.7410
    0.8500, 0.3250, 0.0980
    0.9290, 0.6940, 0.1250
    0.4940, 0.1840, 0.5560
    0.4660, 0.6740, 0.1880
    0.3010, 0.7450, 0.9330
    0.6350, 0.0780, 0.1840];
end