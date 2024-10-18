function [mapped_data,power,fig,FaceColor]=plotPCA(X,Y,FaceColor,d,flag)
% X: m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
% flag: plot ellipses
if nargin<2
    Y = ones(size(X,2),1);
end
if nargin<4
    d = 2;
end
if nargin<5
    flag = 1;
end
U=sort(unique(Y));
% str=cell(1,length(U));

if nargin<3
    FaceColor = defaultColor(length(unique(Y)));
    
end

[mapped_data,~,power]=compute_mapping(X','PCA',d);
mapped_data=mapped_data';

fig = figure;
hold on
if d==3
    for i=1:length(U)
        plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    zlabel(['PC3 (' num2str(round(power(3)*10000)/100) '%)']);
else
    for i=1:length(U)
        plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    if flag == 1
        for i=1:length(U)
            if sum(Y==U(i))>=4
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(i,:))
            end
        end
    end
end
set(gca,'FontSize',14);
% set(gca,'FontSize',14,'FontWeight','Bold');
pbaspect([1 1 1]);
box on
end
