function  FaceColor=ColorSel(s)
if nargin<1
    s=1;
end
switch s
    case 1
        tf=cbrewer('qual', 'Paired',12);
        FaceColor(1,:)=tf(4,:);
        FaceColor(2,:)=tf(10,:);
        FaceColor(3,:)=tf(1,:);
        FaceColor(4,:)=tf(6,:);
        FaceColor(5,:)=tf(5,:);
        FaceColor(6,:)=tf(7,:);
        FaceColor(7,:)=tf(8,:);
    case 2
        tf2=cbrewer('qual', 'Dark2',5);
        tf=cbrewer('qual', 'Paired',12);
        %         FaceColor(1,:)=tf(4,:);
        %         FaceColor(2,:)=tf2(1,:);
        %         FaceColor(3,:)=tf2(2,:);
        %         FaceColor(4,:)=tf2(3,:);
        %         FaceColor(5,:)=tf2(4,:);
        FaceColor(1,:)=tf(4,:);
        FaceColor(2,:)=tf(9,:);
        FaceColor(3,:)=tf(12,:);
        FaceColor(4,:)=tf(8,:);
        FaceColor(5,:)=tf2(4,:);
        FaceColor(6,:)=tf(5,:);
end
end