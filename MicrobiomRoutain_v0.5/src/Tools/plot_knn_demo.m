function plot_knn_demo
clc;close all;clear
[X,Y]=load_Data('spiral');
X=X(1:2,:);
figure,
plot(X(1,:),X(2,:),'o');
grid
D=[X(1,:)' X(2,:)'];
label=Y-1;
[IDX]=plot_knn(D,1,2,'euclidean',label);
end
