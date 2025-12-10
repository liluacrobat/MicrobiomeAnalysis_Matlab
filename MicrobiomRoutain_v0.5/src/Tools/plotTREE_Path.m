function plotTREE_Path(X,F,Y,stree,root,path,path_point)
FaceColor = cbrewer('qual', 'Set2',12);
mapped_data=X;
Y=ones(size(Y))*2;
Y(path_point)=1;
U=sort(unique(Y));
str=cell(1,length(U));
figure,
for i=1:length(U)
    plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:));
    hold on
    str{i}=num2str(U(i));
end
hold on
[m,n]=size(stree);
for i=1:m
    for j=1:n
        if stree(i,j)~=0
            plot3([F(1,i) F(1,j)],[F(2,i) F(2,j)],[F(3,i) F(3,j)],'-k','linewidth',4);
        end
    end
end
plot3(F(1,path),F(2,path),F(3,path),'or','MarkerFaceColor','r');
grid
end