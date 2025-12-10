function StratifyTable(file,level)
clc;close all;clear

if nargin<2
    level=7;
end
file = 'merged_read_table_k2_std.txt';
level = 7;
tb = readtable(file);
reads = table2array(tb(:,2:end));
sample = tb.Properties.VariableNames(2:end);
otu = table2cell(tb(:,1));
x = cell(1,level);
y = cell(1,level);
xu = cell(1,level);
yu = cell(1,level);

flag = zeros(1,level);
p = 1;
tax = [];
for k=1:length(otu)
    t = otu{k};
    s = strsplit(t,'|');
    l = length(s);
    if l==p
        flag(l) = flag(l)+1;
        tax{flag(l)} = s{p};
        x{p}(flag(l),:)=reads(k,:);
    end
end
yall= tax;
xall = x{p};
writeResultTable([file '.k.all.txt'],yall,sample,xall);
for kd = 1:length(yall)
    target = yall{kd};
    switch target
        case 'k__Bacteria'
            level = 7;
        case 'k__Archaea'
            level = 7;
        case 'k__Eukaryota'
            level = 8;
        case 'k__Viruses'
            level = 8;
        otherwise
            level = 0;
    end
    x = cell(1,level);
    y = cell(1,level);
    x{1}(kd,:) = xall(kd,:);
    y{1} = yall;
    xu = cell(1,level);
    yu = cell(1,level);
    for p=2:level
        flag = zeros(1,8);
        yst = zeros(size(x{p-1}));
        tax = [];
        for k=1:length(otu)
            t = otu{k};
            s = strsplit(t,'|');
            l = length(s);
            if strcmp(s{1},target)==1
                flag(l) = flag(l)+1;
                if l==p
                    tax{flag(l)} = s{p};
                    x{p}(flag(l),:)=reads(k,:);
                    fgq = 0;
                    for q = 1:length(y{p-1})
                        if strcmp(y{p-1}{q},s{p-1})==1
                            yst(q,:) = yst(q,:)+reads(k,:);
                            fgq = 1;
                            break
                        end
                    end
                    if fgq ==0
                        disp('Error!!!!!!!!');
                    end
                end
            end
            
        end
        xu{p-1} = x{p-1}-yst;
        y{p} = tax;
    end
    for p=1:level-1
        temp = y{p};
        for i=1:length(temp)
            temp{i} = [temp{i} '_unclassified'];
        end
        yu{p} = temp;
    end
    
    for p=2:level
        xft = x{p};
        yft = y{p}(:);
        nc = 0;
        for q = 1:p-1
            xft = [xft;xu{q}];
            yft = [yft;yu{q}(:)];
            nc = nc+length(sum(xu{q},2)>0);
        end
        tnc(kd,p) = length(y{p});
        tn(kd,p) = nc;
        idx = sum(xft,2)>0;
        xf{p} = xft(idx,:);
        yf{p} = yft(idx);
        if sum(abs(sum(xf{p},1)-sum(x{1},1)))~=0
            disp('Error!!!!!!!!');
        end
        
        writeResultTable(['k2_std/' file '.' target '.L' num2str(p) '.txt'],yf{p},sample,xf{p});
    end
end
for k=1:8
    uL{1,k} = ['Level' num2str(k)];
end
writeResultTable(['k2_std/' file '.all.n_unclassified.txt'],yall,uL,tn);
writeResultTable(['k2_std/' file '.all.n_classified.txt'],yall,uL,tnc);
end