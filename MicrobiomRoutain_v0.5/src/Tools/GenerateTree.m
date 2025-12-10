function [x,y,C]=GenerateTree(Para)
% Generate tree like 2-D data
% Para.cluster: cluster number
% Para.Lmin: min length
% Para.Lmax: max length
% Para.n: number of data in each cluster
% x:sample
% y:label
% C:cluster information
cluster=Para.cluster; 
Lmin=Para.Lmin; 
Lmax=Para.Lmax;
% overLap=Para.overLap;
n=Para.n; 
C=[];
x=rand(2,1)/2;
ang=rand(1)*pi/4;
x=gen_bran(x,n,Lmax,Lmin,ang,pi/2);
as=pi/2+ang;
C=[C x(:,round(n/2))];
for i=1:cluster
    ang=pi/2-(rand(2,1))*pi/8;
%     ang(2)=pi/2-ang(2);
    x1=gen_bran(x(:,end),n,Lmax,Lmin,-ang(1),as+pi/2);
    x2=gen_bran(x(:,end),n,Lmax,Lmin,ang(2),as-pi/2);
    if rand(1)>0.5
        x=[x x1 x2];
        as=as-pi/2+ang(2);
        C=[C x1(:,round(n/2))];
        C=[C x2(:,round(n/2))];
    else
        x=[x x2 x1];
        as=as+pi/2-ang(1);
        C=[C x2(:,round(n/2))];
        C=[C x1(:,round(n/2))];
    end
    
end
Num=size(x,2)/n;
y=zeros(size(x,1),1);
for i=1:Num
    y((i-1)*n+1:i*n,1)=ones(n,1)*i;
end
end
function x=gen_bran(x,n,Lmax,Lmin,ang,ds)
        L=rand(1)*(Lmax-Lmin)+Lmin;
        dang=rand(1,n);
        dang=(dang/sum(dang))*ang;
        dl=rand(1,n);
        dl=dl/sum(dl)*L;
    for j=1:n
        ds=ds+dang(j);
        x=[x x(:,end)+[cos(ds)*dl(j);sin(ds)*dl(j)]];
    end
    x=x(:,2:end);
end