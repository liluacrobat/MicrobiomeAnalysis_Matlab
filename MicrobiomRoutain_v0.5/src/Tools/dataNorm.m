function Y = dataNorm(X,method)
%% normalize the data so that each feature is comparable
switch method
    case '01'
        % Normalize to [0, 1]
        [MIN,~] = min(X,[],2);
        [MAX,~] = max(X,[],2);
        STA = MAX-MIN;
        Y = X-MIN*ones(1,size(X,2));
        Y(STA~=0,:) = Y(STA~=0,:)./(STA(STA~=0,1)*ones(1,size(Y,2)));
        Y(STA==0,:)=0;
    case 'std'
        % Normalize to Z-score
        M = mean(X,2);
        STD=std(X,[],2);
        Y = (X-M*ones(1,size(X,2)))./(STD*ones(1,size(X,2)));
        Y(STD==0,:) = 0;
end
end