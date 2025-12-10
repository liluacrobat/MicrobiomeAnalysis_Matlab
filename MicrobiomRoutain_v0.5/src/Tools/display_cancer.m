function mapped_data=display_cancer(X,Y,Name,str)
FaceColor =  distinguishable_colors(max(Y));
if size(X,1)>2
    type=3;
else
    type=2;
end
switch type
    case 3
        [mapped_data,~]=compute_mapping(X','PCA',3);
        mapped_data=mapped_data';
        figure,
        for i=1:max(Y)
            plot3(mapped_data(1,Y==i),mapped_data(2,Y==i),mapped_data(3,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
            colormap(colorcube);
            hold on
        end
    case 2
        [mapped_data,~]=compute_mapping(X','PCA',2);
        mapped_data=mapped_data';
        figure,
        for i=1:max(Y)
            plot(mapped_data(1,Y==i),mapped_data(2,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
            colormap(colorcube);
            hold on
        end
end
grid
if nargin<3
    legend('Basal','HER2','LA','LB','Normal-like','Normal');
else
    switch Name
        case 0
            legend(str);
        case 1
            legend('Basal');
        case 2
            legend('Basal','Normal');
        case 4
            legend('1','2','3','4','Normal');
        case 5
            legend('Basal','HER2','LA','LB','Normal');
        case 6
            legend('Normal','6','7','8','9','10');
        case 9
            legend('1','2','3','4','HER2','LA','LB','Normal-like','Normal');
        otherwise
            legend('1','2','3','4','5','6','7','8','9','10');
    end
end
end