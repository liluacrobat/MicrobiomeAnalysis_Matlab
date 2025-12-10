function h = plotTone(X,Y,mm)
% X: principal coordinate m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
if nargin<2||isempty(Y)
    Y=ones(size(X,2),1);
end

sel = ~isnan(mm);
X = X(:,sel);
Y=Y(sel);
mm = mm(sel);
mm = mm(:)';
tt = data_norm(mm);
tt = round(tt*99)+1;
Y = tt(:);
ColorG = (cbrewer('seq', 'Blues',100));
ColorG = flip(cbrewer('div', 'RdBu',100));
ColorG(ColorG>1)=1;
ColorG(ColorG<0)=0;
U=sort(unique(Y));
str=cell(1,length(U));

FaceColor=ColorG;

D=size(X,1);
if D>3
    X = X(1:3,:);
    D=3;
end
mapped_data=X;
figure,
switch D
    case 2
        if nargin<3
            h = figure;
        end
        for i=1:length(U)
            plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(U(i),:),'MarkerSize',7,'MarkerEdgeColor',FaceColor(U(i),:));
            hold on
            str{i}=num2str(U(i));
        end
        for i=1:length(U)
            if sum(Y==U(i))>=3
                plotEllipse(mapped_data(1:2,Y==U(i)),FaceColor(U(i),:))
            end
        end
        %         grid
    case 3
        if nargin<3
            h = figure;
        end
        for i=1:length(U)
            plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(U(i),:),'MarkerSize',7,'MarkerEdgeColor',FaceColor(U(i),:));
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