function tax = simpleName(x,mod)
tax = cell(size(x));
if nargin<2
    mod='HOMD_GTDB';
end
switch mod
    case 'HOMD_GTDB'
        for i=1:length(x)
            line = x{i};
            s = strsplit(line, ';');
            [head1,tax1] = rmHead(s{end});
            tax{i} = tax1;
            tax{i} = strrep(tax{i},'_',' ');
        end
    case 'HOMD'
        tax = cell(size(x));
        for i=1:length(x)
            line = x{i};
            s = strsplit(line, ';');
            tax3 = s{end-1};
            tax2= s{end};
            if contains(tax2,'_unclassified')
                tax{i} = strcat(strrep(tax2,'_',' '));
            else
                tax{i} = strcat(strrep(tax3,'_',' '),'_',strrep(tax2,'_',' '));
            end
            tax{i} = strrep(tax{i},'_',' ');
        end
    otherwise
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
end