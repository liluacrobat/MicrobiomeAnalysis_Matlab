function betaAnalysis_by_group(fileName,measure,measure_title,groups_y,header,d,para)
if ~isfield(para,'Ellipse')
    para.Ellipse=1;
end
if ~isfield(para,'Ellipse')
    para.confidence=0.95;
end
if ~isfield(para,'confidence')
    para.confidence=0.95;
end
facecolor = defaultColor(length(header)+length(para.group_header));
d = min(d,3);
if para.pca==0
    dis = measure;
    pcoa = f_pcoa(dis,1);
    mapped_data = pcoa.scores(:,1:d)';
    power = pcoa.expl(:,1);
    plotFIGURE(mapped_data,groups_y,facecolor,para);
    set(gca,'FontSize',14);
    xlabel(['PCo1 (' num2str(round(power(1)*10)/10) '%)']);
    ylabel(['PCo2 (' num2str(round(power(2)*10)/10) '%)']);
    if d==3
        ylabel(['PCo3 (' num2str(round(power(2)*10)/10) '%)']);
    end
else
    plotPCA(measure,groups_y,facecolor,d,para);
end
legend([para.group_header]);
set(gca,'FontSize',14);
pbaspect([1 1 1]);
aa=axis;
aa_x = mean(aa(1:2));
aa_y = mean(aa(3:4));
aa_wx = aa(2)-aa(1);
aa_wy = aa(4)-aa(3);

aa = [aa_x-aa_wx*0.55 aa_x+aa_wx*0.55 aa_y-aa_wy*0.55 aa_y+aa_wy*0.55];
axis(aa);
box on
% hold on
% uy=unique(groups_y);
% for i=1:length(unique(groups_y))
%     id = find(groups_y==uy(i))
% plot(mapped_data(1,id)',mapped_data(2,id)','.-k');
% end
keyboard
if para.pca==1
    plotPDF(gcf,strcat(fileName,'_PCA'));
else
    plotPDF(gcf,strcat(fileName,'_PCoA'));
end

end
function h = plotFIGURE(X,Y,FaceColor,para)
% X: principal coordinate m*n data, m rows n samples
% Y: group label, 1,2, 3, 4, 5, 6 ...
% FaceColor: color used for each group
if nargin<2||isempty(Y)
    Y=ones(size(X,2),1);
end


U=sort(unique(Y));
str=cell(1,length(U));
if nargin<3
    FaceColor=cbrewer('qual', 'Set1',9);
end
D=size(X,1);
if D>3
    X = X(1:3,:);
    D=3;
end
mapped_data=X;
U2 = unique(para.group_y);
h = figure;
hold on;

switch D
    case 2

        if para.Ellipse==1
            for i=1:length(U2)
                plot(mapped_data(1,para.group_y==U2(i)),mapped_data(2,para.group_y==U2(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
            end
            for i=1:length(U)
                plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i+length(U2),:),'MarkerSize',8,'MarkerEdgeColor','k');

                str{i}=num2str(U(i));
            end
            for i=1:length(U2)
                if sum(para.group_y==U2(i))>=3
                    plotEllipse(mapped_data(1:2,para.group_y==U2(i)),FaceColor(i,:),para.confidence)
                end
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
pbaspect([1 1 1]);
end
function [mapped_data,power,fig,FaceColor]=plotPCA(X,Y,FaceColor,d,para)
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
U=sort(unique(Y));
str=cell(1,length(U));

if nargin<3
    FaceColor = defaultColor(length(unique(Y)));

end

[mapped_data,~,power]=compute_mapping(X','PCA',d);
mapped_data=mapped_data';

fig = figure;
hold on
U2 = unique(para.group_y);
if d==3
    for i=1:length(U)
        plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    end
    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
    zlabel(['PC3 (' num2str(round(power(3)*10000)/100) '%)']);
else
    if para.Ellipse==1
        for i=1:length(U2)
            plot(mapped_data(1,para.group_y==U2(i)),mapped_data(2,para.group_y==U2(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
        end
        for i=1:length(U)
            plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i+length(U2),:),'MarkerSize',8,'MarkerEdgeColor','k');
        end
        for i=1:length(U2)
            if sum(para.group_y==U2(i))>=3
                plotEllipse(mapped_data(1:2,para.group_y==U2(i)),FaceColor(i,:),para.confidence)
            end
        end
    end

    xlabel(['PC1 (' num2str(round(power(1)*10000)/100) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*10000)/100) '%)']);
end
set(gca,'FontSize',14);
% set(gca,'FontSize',14,'FontWeight','Bold');
pbaspect([1 1 1]);
box on
end
