function plotTREE(Xor,CENor,X,CEN,p,Y)
U=unique(sort(Y));
FaceColor =  distinguishable_colors(length(U));
if nargin<5
    p=1;
end
figure,
Y=Y(1:length(Y)/2);
hold on
[IDX,dis] = knnsearch(CENor',Xor','K',1);
for i=1:length(U)
    temp=find(Y==U(i));
plot3([p*X(1,Y==U(i))+(1-p)*CEN(1,IDX(temp))],[p*X(2,Y==U(i))+(1-p)*CEN(2,IDX(temp))],[p*X(3,Y==U(i))+(1-p)*CEN(3,IDX(temp))],'o','MarkerFaceColor',FaceColor(U(i),:));
end
plot3([CEN(1,:)],[CEN(2,:)],[CEN(3,:)],'ok');
for i=1:length(U)
    name{i}=num2str(U(i));
end
name{end}='tree';
legend(name);
for i=1:length(IDX)
plot3([p*X(1,i)+(1-p)*CEN(1,IDX(i)) CEN(1,IDX(i))],[p*X(2,i)+(1-p)*CEN(2,IDX(i)) CEN(2,IDX(i))],[p*X(3,i)+(1-p)*CEN(3,IDX(i)) CEN(3,IDX(i))],'-','color',FaceColor(Y(i),:));
end
hold off
grid;

end