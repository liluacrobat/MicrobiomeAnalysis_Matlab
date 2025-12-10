function [x,y]=GenTree_data(v,Para)
if nargin<2
    Para.cluster=2; 
    Para.Lmin=0.2; 
    Para.Lmax=1;
    Para.n=200; 
    Para.C=0.05;
    Para.kn=10;
end
if nargin<1
    v=0;
end
% [x1,y1,C1]=GenerateTree(Para);
% [x2,y2,C2]=GenerateTree(Para);
% [x3,y3,C3]=GenerateTree(Para);
load tree_x123 x1 x2 x3 y1 y2 y3
sig=0.1;
sign=0.5;
rng default
x1=x1+randn(size(x1))*sig;
x2=x2+randn(size(x2))*sig;
x3=x3+randn(size(x2))*sig;
tep=x3;
x3=x1;
x1=tep;
tep=y3;
y3=y1;
y1=tep;
r=rand(1,size(x1,2));
[~,id]=sort(r);
x1(:,:)=x1(:,id);
y1=y1(id);
x1=add_cluster(x1,y1,Para.C);
x23=add_cluster([x2],y2,Para.C);
[x23t,y2]=add_noise(x23,y2,Para.n,Para.kn);
x23(size(x23,1)/2+1:end,:)=x23t(size(x23,1)/2+1:end,:);
x=[x1;x23];
% x(5,:)=1+0.001*rand(1,length(x(7,:)));
% x(6,:)=0+0.01*rand(1,length(x(7,:)));
x=data_norm(x);
if v==1
display_cancer(x,y2);
title('Structure Data');
display_cancer(x(1:2,:),y2);
title('Structure 1:Tree');
display_cancer(x(3:4,:),y2);
title('Structure 1:Cluster');
display_cancer(x(1:4,:),y2);
title('Structure 1');
temp=x(5:end,:);
l=size(temp,1)/2;
display_cancer(temp(1:l,:),y2);
title('Structure 2:Tree');
display_cancer(temp(l+1:end,:),y2);
title('Structure 2:Cluster');
display_cancer(temp,y2);
title('Structure 2');
end
for i=1:1000
    x=[x;randn(1,size(x1,2))*sign];
end
r=rand(1,size(x,2));
[~,id]=sort(r);
x(13:14,id(10:end))=x(13:14,id(10:end))/100;
[x,y2]=randSel(x,y2,1);
y=y2;
% save toy_data_tree_refine_lasso2 x y1 y2 
end
function [x,y]=add_noise(x,y,n,kn)
for i=n:2*n:length(y)-n-kn
    temp=[y(i-kn:i+kn);y(i+n:i+n+kn)];
    tx=[x(:,i-kn:i+kn) x(:,i+n:i+n+kn)];
    Perm=randperm(length(temp));
    y([i-kn:i+kn i+n:i+n+kn])=temp(Perm);
    x(:,[i-kn:i+kn i+n:i+n+kn])=tx(:,Perm);
end
end
function x=add_cluster(x,y2,sig)
x2=x;
for i=1:max(y2)
    temp=randn(size(x(:,y2==i)))*sig;
    C=mean(x(:,y2==i),2);
    x2(:,y2==i)=temp+C*ones(1,size(temp,2));
end
x=[x;x2];
end
function [xt,yt]=randSel(x,y,p)
l=round(length(y)*p);
rng default
t=randn(1,length(y));
[~,id]=sort(t);
yt=y(id(1:l));
xt=x(:,id(1:l));
end