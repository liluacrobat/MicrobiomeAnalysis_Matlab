function [y,y2] = refineOTUName(x)
for i=1:length(x)
    temp = x{i};
    str_ls = strsplit(temp, ';');
    n = length(str_ls);
    if n==1
        ns = str_ls{1};
        ns2 = ns;
    else
        otu_ls = [];
        for j=1:n
            
            t = strtrim(str_ls{j});
            if strcmp(t(3),'_')==1
                otu_ls{j} = t(4:end);
            end
        end
        for j=n:-1:1
            if strcmp(otu_ls{j},'unclassified')~=1
                if n==7
                    ns = [otu_ls{6} ' ' otu_ls{7}];
                else
                    ns = [otu_ls{j} ' unclassified'];
                end
            end
        end
        ns2 = otu_ls{2};
    end
    ns = strrep(ns,'_',' ');
    y{i} = ns;
    y2{i} = ns2;
end