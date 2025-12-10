function plotResult(OTU_ID,Rel,Type,M,Q,alpha,Mthr)
U=unique(Type);
ids=M>Mthr&Q<alpha;
P_Rel=Rel(ids==1,Type==U(2));
N_Rel=Rel(ids==1,Type==U(1));
Ab_ID=OTU_ID(ids==1);
P_M=mean(P_Rel,2);
N_M=mean(N_Rel,2);
diff=P_M-N_M;

[diff,idx]=sort(diff);
P_Rel=P_Rel(idx,:);
N_Rel=N_Rel(idx,:);
Ab_ID=Ab_ID(idx);
% 
% P_id=Ab_ID(diff>0);
% P_Rel=P_Rel(diff>0,:);
% N_id=Ab_ID(diff<0);
% N_Rel=N_Rel(diff<0,:);
figure,
hold on
for i=1:length(Ab_ID)
    boxplot(N_Rel(i,:)',[],'Positions',i-0.5,'Colors','g','Widths',0.3);
end
boxplot(N_Rel','Positions',[1:length(Ab_ID)]-0.5,'Colors','g','Widths',0.3);
% boxplot(N_Rel','Label',Ab_ID,'Colors','g');
hold on
boxplot(P_Rel','positions',[1:length(Ab_ID)]+0.5,'Colors','r');


xticks([1:length(Ab_ID)]);
set(gca,'XTickLabel', Ab_ID);
xtickangle(-30)

end