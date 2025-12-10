function [DATA,STD]=data_norm2(DATA)
%% normalize the data so that each feature is comparable
% Normalize to mean 0 std 1
M=mean(DATA,2);
Q = quantile(DATA,[0.05 0.95],2);
for i=1:size(DATA,1)
    tmp = DATA(i,:);
    sel = tmp>Q(1) & tmp<Q(2);
    STD(i,1) = std(tmp);
end
DATA=(DATA-M*ones(1,size(DATA,2)))./(STD*ones(1,size(DATA,2)));
DATA(STD==0,:)=0;
end
