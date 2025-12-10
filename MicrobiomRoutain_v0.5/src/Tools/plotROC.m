function [AUC, roc_x, roc_y]= plotROC(Y,Score)
UY = unique(Y);
N = length(UY);
if size(Score,2)==1
    Score(:,2) = Score;
end
Facecolor =defaultColor;
fig = figure;
hold on;
if N==2
    [roc_x, roc_y,T,AUC,OPTROCPT] = perfcurve(Y,Score(:,2),max(UY));
    plot(roc_x, roc_y,'linewidth',2,'color','k');
else
    for i=1:N
        [roc_x, roc_y,T,AUCt,OPTROCPT] = perfcurve(Y,Score(:,i),UY(i));
        AUC_seq(i) = AUCt;
        plot(roc_x, roc_y,'linewidth',2,'color',Facecolor(mod(i-1,7)+1,:));
    end
    
    AUC = mean(AUC_seq);
end
xlabel('1-Specificity');
ylabel('Sensitivity');
box on
set(gca,'FontSize',14);
pbaspect([1 1 1])
end