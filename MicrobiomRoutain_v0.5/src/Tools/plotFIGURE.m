function h = plotFIGURE(X,Y,FaceColor)
% X: principal coordinate m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
if nargin<2||isempty(Y)
    Y=ones(size(X,2),1);
end


U=sort(unique(Y));
str=cell(1,length(U));
if nargin<3
    % FaceColor=cbrewer('qual', 'Set1',9);
    FaceColor=defaultColor(length(unique(Y)));
end
D=size(X,1);
if D>3
    X = X(1:3,:);
    D=3;
end
mapped_data=X;
switch D
    case 2
        if nargin<3
            h = figure;
        end
        for i=1:length(U)
            plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
            hold on
            str{i}=num2str(U(i));
        end
        for i=1:length(U)
            if sum(Y==U(i))>=3
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(i,:))
            end
        end
        %         grid
    case 3
        if nargin<3
            h = figure;
        end
        for i=1:length(U)
            plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
            hold on
            str{i}=num2str(U(i));
        end
        
        grid
end
% boldify
set(gca,'FontSize',14);
% set(gca,'FontSize',14,'FontWeight','Bold');
% pbaspect([1 1 1]);
end