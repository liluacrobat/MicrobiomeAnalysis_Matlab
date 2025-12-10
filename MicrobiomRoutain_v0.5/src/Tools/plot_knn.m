function [IDX]=plot_knn(D,k,dim,dis,label)
 % plot the graph of knn
 if nargin<2
     k=1;
 end
 if nargin<3
     dim=3;
 end
 if nargin<4
     dis='euclidean';
 end
 if nargin<5
     label=zeros(size(D,1),1);
 end
IDX=knnsearch(D,D,'K',k+1,'Distance',dis);
[PCA_D,~]=compute_mapping(D,'PCA',dim);
colormap(colorcube);
figure,
i=0;
index=find(label==i);
MARK='-or';
plot_pair(PCA_D,IDX,index',MARK,dim);
grid
hold on
i=1;
index=find(label==i);
MARK='-ob';
plot_pair(PCA_D,IDX,index,MARK,dim)
hold on
i=2;
index=find(label==i);
MARK='-oc';
plot_pair(PCA_D,IDX,index,MARK,dim)
hold on
i=3;
index=find(label==i);
MARK='-ok';
plot_pair(PCA_D,IDX,index,MARK,dim)
hold on

% if dim==3
%     plot3(PCA_D(:,1)',PCA_D(:,2)', PCA_D(:,3)','ob');
% else
%     plot(PCA_D(:,1)',PCA_D(:,2)','ob');
% end
grid
colormap(colorcube);
hold off
end