function FaceColor=ColorList(name)
switch name
    case 'ibd_pca_be'
        tf=cbrewer('qual', 'Paired',12);
        FaceColor(1,:)=tf(4,:);
        FaceColor(2,:)=tf(10,:);
        FaceColor(3,:)=tf(1,:);
        FaceColor(4,:)=tf(6,:);
        tf=FaceColor;
        tf2 = cbrewer('qual', 'Set2',8);
        FaceColor=tf2;
        FaceColor(2:4,:)=tf2(2:4,:);
        FaceColor(1,:) =  tf(1,:);
        
    case 'ibd_pca_be2'
        tf=cbrewer('qual', 'Paired',12);
        FaceColor(1,:)=tf(4,:);
        FaceColor(2,:)=[219 109 0]/255;
        FaceColor(3,:)=[182 219 255]/255;%[182 109 255]/255;
        FaceColor(4,:)=[255 182 219]/255;
       
    case 'ibd_ddr_be'
        tf=cbrewer('qual', 'Set3',12);
        tf2=cbrewer('qual', 'Set1',9);
        FaceColor = [tf2(3,:);
            tf(4,:);
            tf(10,:);
            tf(9,:)];
        
        FaceColor = [tf2(3,:);
            tf(4,:);
            tf(10,:);
            tf(9,:)];
        tf2 = cbrewer('qual', 'Set2',8);
        FaceColor=tf2;
        FaceColor(3:5,:)=tf2(2:4,:);
        FaceColor(1,:) =  tf2(1,:);
        FaceColor(2,:) =  tf2(5,:);
        tf3 = cbrewer('qual', 'Set1',9);
        tf4 = cbrewer('qual', 'Set3',12);
        FaceColor(1,:) =  tf3(3,:);
        FaceColor = [FaceColor(1,:); FaceColor(3:end,:)];
        FaceColor(4,:) = tf4(9,:);
    case  'ibd_ddr_dia'
        tf=cbrewer('qual', 'Paired',12);
        FaceColor(3,:)=tf(10,:);
        FaceColor(4,:)=tf(1,:);
        FaceColor(5,:)=tf(6,:);
        tf2 = cbrewer('qual', 'Set2',8);
        FaceColor(1,:) =  tf2(1,:);
        FaceColor(2,:) =  tf2(5,:);
    case  'ibd_pca_dia'
        tf=cbrewer('qual', 'Paired',12);
        FaceColor(2,:)=tf(10,:);
        FaceColor(3,:)=tf(1,:);
        FaceColor(4,:)=tf(6,:);
        FaceColor(1,:) =  tf(4,:);
        case  'ibd_pca_dia2'
        tf=cbrewer('qual', 'Paired',12);
        FaceColor(2,:)=tf(10,:);
        FaceColor(3,:)=tf(12,:);
        FaceColor(4,:)=tf(1,:);
        FaceColor(5,:)=tf(5,:);
        FaceColor(6,:)=tf(6,:);
        FaceColor(1,:) =  tf(4,:);
    case 'ibd_cluster'
        tf2 = cbrewer('qual', 'Set2',8);
        FaceColor(1,:) =  tf2(1,:);
        FaceColor(2,:) =  tf2(5,:);
        tf = cbrewer('qual', 'Set3',12);
        FaceColor(3:6,:) =  tf(3:6,:);
        FaceColor(7:9,:) =  tf(8:10,:);
        FaceColor(9,:) =  tf(11,:);
        FaceColor(1,:) =  tf(1,:);
        FaceColor(2,:) =  tf(7,:);
        FaceColor(3,:) =  tf(3,:);
        FaceColor(4,:) =  tf(6,:);
        FaceColor(5,:) =  tf(4,:);
%                 FaceColor = FaceColor(2:end,:);
%                 FaceColor = [FaceColor(1,:);FaceColor(4:6,:);FaceColor(8,:)];
%                 FaceColor = tf2;
    otherwise
        FaceColor =  cbrewer('qual', 'Set1',length(unique(Y)));
        
end
end
