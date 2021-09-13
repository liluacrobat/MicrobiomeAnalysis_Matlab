function tax = simpleName(x)
tax = cell(size(x));
for i=1:length(x)
    line = x{i};
    s = strsplit(line, ';');
    [head1,tax1] = rmHead(s{end});
    
    if strcmp(head1,'s')==1
        [head2,tax2] = rmHead(s{end-1});
        tax{i} = strcat(strrep(tax2,'_',' '),'_',strrep(tax1,'_',' '));
    else
        tax{i} = strrep(tax1,'_',' ');
    end
    tax{i} = strrep(tax{i},'_',' ');
end
end