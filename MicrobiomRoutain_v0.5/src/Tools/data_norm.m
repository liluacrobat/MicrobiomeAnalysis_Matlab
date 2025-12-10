function DATA=data_norm(DATA)
%% normalize the data so that each feature is comparable
% Normalize to [0, 1]
[MIN,~] = min(DATA,[],2);
[MAX,~] = max(DATA,[],2);
STA=MAX-MIN;
DATA=DATA-MIN*ones(1,size(DATA,2));
DATA(STA~=0,:)=DATA(STA~=0,:)./(STA(STA~=0,1)*ones(1,size(DATA,2)));
DATA(STA==0,:)=0;
end