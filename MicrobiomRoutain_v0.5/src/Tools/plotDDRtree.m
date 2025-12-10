function plotDDRtree(Z,tr,stree,Y,Legend,flag,FaceColor)
if nargin<6
    flag=1;
end
figure,
mapped_data=Z;
if nargin<7
    FaceColor=ColorSel(1);
    switch length(unique(Y))
        case 6
            temp=FaceColor;
            FaceColor(4,:)=temp(6,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 3
            temp=FaceColor;
            FaceColor(3,:)=temp(4,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 2
            temp=FaceColor;
            FaceColor(2,:)=temp(4,:);
            FaceColor(5:6,:)=temp(4:5,:);
        case 4
            tf=cbrewer('qual', 'Paired',12);
            
            FaceColor(1,:)=tf(4,:);
            FaceColor(2,:)=tf(10,:);
            FaceColor(3,:)=tf(1,:);
            FaceColor(4,:)=tf(6,:);
            
            
        otherwise
            tf=cbrewer('qual', 'Paired',12);
            FaceColor(1,:)=tf(4,:);
            FaceColor(2,:)=tf(10,:);
            FaceColor(3,:)=tf(1,:);
            FaceColor(4,:)=tf(6,:);
            tf=FaceColor;
            tf2 = cbrewer('qual', 'Set2',length(unique(Y)));
            FaceColor =  cbrewer('qual', 'Set2',length(unique(Y)));
            FaceColor(3:5,:)=tf(2:4,:);
    end
end
%plot smaples


grid;
U=unique(Y);
if size(tr,1)>=3
    for i=1:length(U)
        if ~isempty(find(Y==U(i), 1))
            plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',10,'MarkerEdgeColor','k');
            colormap(colorcube);
            hold on
        end
    end
    PCA_D=tr(1:3,:);
    hold on
    [m,n]=size(stree);
    % plot curve
    if flag~=0
        for i=1:m
            for j=1:n
                if stree(i,j)~=0
                    plot3([PCA_D(1,i) PCA_D(1,j)],[PCA_D(2,i) PCA_D(2,j)],[PCA_D(3,i) PCA_D(3,j)],'-k','linewidth',4);
                end
            end
        end
    end
    % plot3(PCA_D(1,root_index),PCA_D(2,root_index),PCA_D(3,root_index),'or','MarkerFaceColor','r','MarkerSize',20);
    IDX = knnsearch(tr',Z','K',1);
    m=size(Z,2);
    Para.knn=1;
    for i=1:m
        for j=1:Para.knn
            %         i
            %         j
            plot3([Z(1,i) tr(1,IDX(i,j))],[Z(2,i) tr(2,IDX(i,j))],[Z(3,i) tr(3,IDX(i,j))],'--','linewidth',2,'color',FaceColor(find(Y(i)==U),:));
        end
    end
    xlabel('DDR1')
    ylabel('DDR2')
    zlabel('DDR3')
    legend(Legend);
else
    for i=1:length(U)
        if ~isempty(find(Y==U(i), 1))
            plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor',FaceColor(i,:));
            colormap(colorcube);
            hold on
        end
    end
    PCA_D=tr(1:2,:);
    hold on
    [m,n]=size(stree);
    % plot curve
    for i=1:m
        for j=1:n
            if stree(i,j)~=0
                plot([PCA_D(1,i) PCA_D(1,j)],[PCA_D(2,i) PCA_D(2,j)],'-k','linewidth',3);
            end
        end
    end
    % plot3(PCA_D(1,root_index),PCA_D(2,root_index),PCA_D(3,root_index),'or','MarkerFaceColor','r','MarkerSize',20);
    IDX = knnsearch(tr',Z','K',1);
    m=size(Z,2);
    Para.knn=1;
    for i=1:m
        for j=1:Para.knn
            plot([Z(1,i) tr(1,IDX(i,j))],[Z(2,i) tr(2,IDX(i,j))],'--','linewidth',1,'color',FaceColor(find(Y(i)==U),:));
        end
    end
    xlabel('DDR1')
    ylabel('DDR2')
    legend(Legend);
end
end