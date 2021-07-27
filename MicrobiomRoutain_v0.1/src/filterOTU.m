function sel = filterOTU(rel,para)

if isfield(para,'thr_fre')
    para.thr_fre = 0;
end
if isfield(para,'thr_rel')
    para.thr_frerel = 0;
end
y = para.y;
u = unique(y);
n = length(u);
d = size(rel,1);
fre = zeros(d,n);
averel = zeros(d,n);

for i=1:n
    temp = rel(:,y==u(i));
    fre(:,i) = sum(temp>0,2)/sum(y==u(i));
    averel(:,i) = mean(temp,2);
end

max_fre = max(fre,[],2);
max_averel = max(averel,[],2);
sel = max_fre>para.thr_fre & max_averel>para.thr_rel;
end