function [AUC, roc_x, roc_y, AUC_c]= CalAUC(Y,Score)
UY = unique(Y);
N = length(UY);
if size(Score,2)==1
    Score(:,2) = Score;
end
if N==2
    [roc_x, roc_y,T,AUC,OPTROCPT] = perfcurve(Y,Score(:,2),max(UY));
    AUC_c = AUC;
else
    for i=1:N
        [roc_x, roc_y,T,AUCt,OPTROCPT] = perfcurve(Y,Score(:,i),UY(i));
        AUC_seq(i) = AUCt;
    end
    AUC_c = AUC_seq;
    AUC = mean(AUC_seq);
end
end