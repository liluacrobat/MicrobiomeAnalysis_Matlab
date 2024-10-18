function dis = loadBeta(fid,sid)
tbl = readtable(fid,'delimiter','\t','FileType','delimitedtext');
s = table2array(tbl(:,1));
[idx1,~] = AlignID(sid,s);
head = sid(idx1==0);
s_h = [head;s(:)];
if length(s_h)~=length(sid)
    error('Sample list mismatch');
end
[idx2,~] = AlignID(sid,s_h);
n = length(sid);
dis_l = zeros(n);
[n1,n2] = size(tbl);
for i=2:n2
    temp = table2array(tbl(:,i));
    if iscell(temp)
        tmp1 = zeros(size(temp));
        for j=1:length(tmp1)
            pdis = str2num(temp{j});
            if isempty(pdis)
                tmp1(j) = nan;
            else
                tmp1(j) = pdis;
            end
        end
    else
        tmp1 = temp;
    end
    dis_l(2:end,i-1) = tmp1;
end
% dis_l(2:end,1:end-1) = table2array(tbl(:,2:end));
dis_l(isnan(dis_l)) = 0;
dis_full = dis_l+dis_l';
dis = dis_full(idx2,idx2);
end