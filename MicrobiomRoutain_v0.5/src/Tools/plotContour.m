function [mapped_data,power]=plotContour(X,Y,FaceColor)
if nargin<2
    Y=ones(size(X,2),1);
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


[mapped_data,~,power]=compute_mapping(X','PCA',size(X,1));
mapped_data=mapped_data';
for i=1:length(U)
    temp=mapped_data(1:3,Y==U(i));
    %     temp=X(:,Y==U(i));
    Center=mean(temp,2);
    dis=pdist2(temp',Center');
    [~,idx]=sort(dis,'descend');
    sel_idx{i}=dis<dis(idx(round(size(temp,2)*0.1)));
end
figure,
for i=1:length(U)
    temp=mapped_data(:,Y==U(i));
    plot3(temp(1,sel_idx{i}==1),temp(2,sel_idx{i}==1),temp(3,sel_idx{i}==1),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
    hold on
    str{i}=num2str(U(i));
end

xlabel(['PC1 (' num2str(round(power(1)*1000)/10) '%)']);
ylabel(['PC2 (' num2str(round(power(2)*1000)/10) '%)']);
zlabel(['PC3 (' num2str(round(power(3)*1000)/10) '%)']);

grid
boldify_line;
for i=1:length(U)
    temp=mapped_data(:,Y==U(i));
    temp=temp(:,sel_idx{i}==1);
    plotEllipsoid(temp(1:3,:),FaceColor(i,:))
end
legend(str)
end
