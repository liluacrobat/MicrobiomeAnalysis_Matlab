function [mapped_data,power]=plotPCA3D(X,Y)
if nargin<2
    Y=ones(size(X,2),1);
end
U=sort(unique(Y));
str=cell(1,length(U));
FaceColor =  distinguishable_colors(length(U));
n=size(X,1);
if n<3
    X=[X;zeros(3-n,size(X,2)) ];
end
    
[mapped_data,~,power]=compute_mapping(X','PCA',size(X,1));
mapped_data=mapped_data';
figure,
for i=1:length(U)
    plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8);
    colormap(colorcube);
    hold on
    str{i}=num2str(U(i));
end
xlabel(['PC1 (' num2str(round(power(1)*1000)/10) '%)']);
ylabel(['PC2 (' num2str(round(power(2)*1000)/10) '%)']);
zlabel(['PC3 (' num2str(round(power(3)*1000)/10) '%)']);
    legend(str)
grid
end