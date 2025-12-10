function mapped_data=plotPCA2D(X,Y)
if nargin<2
    Y=ones(size(X,2),1);
end
U=sort(unique(Y));
str=cell(1,length(U));
FaceColor =  distinguishable_colors(length(U));
[mapped_data,~]=compute_mapping(X','PCA',2);
mapped_data=mapped_data';
figure,
for i=1:length(U)
    plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:));
    colormap(colorcube);
    hold on
    str{i}=num2str(U(i));
end
    legend(str)
grid
end