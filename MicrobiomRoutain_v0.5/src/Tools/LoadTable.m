function X = LoadTable(filename)
fid = fopen(filename,'r');
tline = fgetl(fid);
sline = strsplit(tline,'\t');
X = cell(1,length(sline));
X(1,:) = sline;
n = 1;
while ~feof(fid)
    n = n+1;
    tline = fgetl(fid);
    sline = strsplit(tline,'\t');
    X(n,1:length(sline)) = sline;
end
end