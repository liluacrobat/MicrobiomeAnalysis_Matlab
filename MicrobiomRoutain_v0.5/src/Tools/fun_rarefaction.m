function X_s = fun_rarefaction(X,depth)
%% rarefact OTU table into a specified depth
% X: original OTU table
% depth: target rarefaction depth
% X_s: rarefacted OTU table
rng default
% cumulated depth
X = round(X);
[N,M]=size(X);
C = round(cumsum(X,1));
% expand each otu table into independent observations
Xfull = cell(1,M);
for i = 1:size(X,2)
    temp = zeros(C(end,i),1);
    temp(1:C(1,i)) = 1;
    for j = 2:N
        temp(C(j-1,i)+1:C(j,i)) = j;
    end
    Xfull{i} = temp;
end

X_s = zeros(size(X));
for i = 1:size(X,2)
    total = max(C(end,i),depth);
    X_s(:,i) = fun_pick(total,depth,Xfull{i},N);
end

end
function X_s = fun_pick(total,dp,expandTable,N)
rng default
X_s = zeros(N,1);
sel = randperm(total,dp); % generate num random values from permutation
idx = zeros(total,1);
idx(sel) = 1; % index of the picked
temp = expandTable;
temp(idx==0) = 0; % Unpicked are assigned as 0
u = setdiff(unique(temp),0); % Non-zero OTU idex
for k = 1:length(u)
    X_s(u(k),1) = sum(temp==u(k));
end
end