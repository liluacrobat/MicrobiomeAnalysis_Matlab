function [m1,m2,LogFC,P,Q] = CalPairGroup(X,Y,method)
U = unique(Y);
LogFC = zeros(size(Y(:)));
P = zeros(size(Y(:)));
Q = zeros(size(Y(:)));
if nargin<3
    method = 'rank';
end
if length(U)==2
    m1 = mean(X(:,Y==U(1)),2);
    m2 = mean(X(:,Y==U(2)),2);
    d = size(X,1);
    for i=1:d
        temp = X(i,:);
        LogFC(i,1) = CalFC(temp,Y);
        switch method
            case 't-test'
                [~,P(i,1)] = ttest2(temp(Y==U(1)),temp(Y==U(2)));
            case 'rank'
                [P(i,1),~] = ranksum(temp(Y==U(1)),temp(Y==U(2)));
        end
    end
    Q = fwer_bonf(P, 0.05, 'off');
    if min(X(:))>=0
        for i=1:d
            temp = X(i,:);
        end
    end
    
else
    disp('Error...')
    
end

end