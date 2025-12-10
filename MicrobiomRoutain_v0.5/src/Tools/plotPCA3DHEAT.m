function plotPCA3DHEAT(X,Y,N,Label,logflag)
if nargin<3
    N=10;
end
if nargin<5
    logflag=1;
end
EdgeColor =  distinguishable_colors(length(unique(Label)));
X=X(:,~isnan(Y));
Y=Y(~isnan(Y));
switch logflag
    case 1
U=logspace(log10(min(Y)), log10(max(Y)),N+1);
sU=logspace(log10(min(Y)), log10(max(Y)),11);
    case 0
U=linspace((min(Y)), (max(Y)),N+1);
sU=linspace((min(Y)), (max(Y)),11);
end
FaceColor =  jet(N);
n=size(X,1);
if n<3
    X=[X;zeros(3-n,size(X,2)) ];
end
    
[mapped_data,~,power]=compute_mapping(X','PCA',3);
mapped_data=mapped_data';
figure,
if nargin>3
UU=unique(Label);
for i=1:length(UU)
    plot3(mapped_data(1,Label==UU(i)),mapped_data(2,Label==UU(i)),mapped_data(3,Label==UU(i)),'o','MarkerSize',12,'MarkerFaceColor',EdgeColor(i,:),'LineWidth',1);
    hold on
end
end
for i=1:length(U)-1
    plot3(mapped_data(1,Y>=U(i)&Y<=U(i+1)),mapped_data(2,Y>=U(i)&Y<=U(i+1)),mapped_data(3,Y>=U(i)&Y<=U(i+1)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','none');
    hold on
end

colormap(FaceColor);
colorbar('TickLabels',round(sU));
xlabel(['PC1 (' num2str(round(power(1)*1000)/10) '%)']);
ylabel(['PC2 (' num2str(round(power(2)*1000)/10) '%)']);
zlabel(['PC3 (' num2str(round(power(3)*1000)/10) '%)']);
grid
end