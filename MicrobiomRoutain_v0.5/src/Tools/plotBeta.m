function [mapped_data,power,fig,FaceColor]=plotBeta(X,Y,FaceColor,flag)
% X: beta dievristy coordinate m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
% flag: plot ellipses
if nargin<2
    Y = ones(size(X,2),1);
end
if nargin<4
    flag = 1;
end
U=sort(unique(Y));
str=cell(1,length(U));
% FaceColor =  distinguishable_colors(length(U));
if nargin<3
    
    FaceColor=ColorSel(1);
    switch length(unique(Y))
        case 6
            temp=FaceColor;
            FaceColor(4,:)=temp(6,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 4
            tf=cbrewer('qual', 'Paired',12);
            FaceColor(1,:)=tf(4,:);
            FaceColor(2,:)=tf(10,:);
            FaceColor(3,:)=tf(1,:);
            FaceColor(4,:)=tf(6,:);
        case 5
            tf=cbrewer('qual', 'Paired',12);
            FaceColor(1,:)=tf(4,:);
            FaceColor(2,:)=tf(10,:);
            FaceColor(3,:)=tf(1,:);
            FaceColor(4,:)=tf(6,:);
            FaceColor(5,:)=tf(8,:);
        case 3
            temp=FaceColor;
            FaceColor(3,:)=temp(4,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 2
            temp=FaceColor;
            FaceColor(2,:)=temp(4,:);
            FaceColor(5:6,:)=temp(4:5,:);
        otherwise
            FaceColor =  cbrewer('qual', 'Set1',length(unique(Y)));
    end
end
n=size(X,1);
if n<3
    X=[X;zeros(3-n,size(X,2)) ];
end
[mapped_data,~,power]=compute_mapping(X','PCA',3);
% [mapped_data,~,power]=compute_mapping(X','PCA',size(X,1));
mapped_data=mapped_data';

fig = figure;
hold on
for i=1:length(U)
    plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
end
xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
zlabel(['PC3 (' num2str(round(power(3)*10000)/100) '%)']);
%     legend(str)
grid
if flag == 1
    for i=1:length(U)
        if sum(Y==U(i))>4
            plotEllipse(mapped_data(1:3,Y==U(i)),FaceColor(i,:))
        end
    end
end

set(gca,'FontSize',14);
% set(gca,'FontSize',14,'FontWeight','Bold');
pbaspect([1 1 1]);
box on

end
