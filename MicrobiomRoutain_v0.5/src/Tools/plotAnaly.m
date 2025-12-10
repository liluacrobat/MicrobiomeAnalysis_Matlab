function plotAnaly(rX,Y,w,k)
[~,id]=sort(w,'descend');
if nargin<4
    X=rX;
else
    X=rX(id(1:k),:);
end
tX=[];
tY=[];
mark=[1 2 3 4 5 6];
for k=1:length(mark)
    i=mark(k);
    tX=[tX,X(:,Y==i)];
    tY=[tY;Y(Y==i)];
end
figure,semilogx(w);
display_cancer(X,Y);
show_stru(tX,tY,'DDR');
% show_stru(tX,tY,'PPT');
% cluster_analy(tX,tY,'all');
% heatmap_analy(tX);
% plotSpanningTree(tX,tY);
% display_cancer(tX,tY);
% cluster_analy(tX,tY,'colum');
% cluster_analy(tX(:,tY==1),tY(tY==1),'all');
% cluster_analy(tX(:,tY==2),tY(tY==2),'all');
% plot_knn(X,Y,3);
% show_stru(tX,tY,'tSNE');
% show_stru(tX,tY,'LLE');
% display_cancer(X,Y);
% show_stru(X,Y,'DDR');
% show_stru(X,Y,'LDA');
% plot_knn(X,Y,1);
end
function [IDX]=plot_knn(X,Y,k)
% plot the graph of knn
N=size(X,2);
dis='cityblock';
IDX=knnsearch(X',X','K',k+1,'Distance',dis);
[PCA_D,~]=compute_mapping(X','PCA',3);
PCA_D=PCA_D';
FaceColor =  distinguishable_colors(max(Y));
figure,
for i=1:max(Y)
    plot3(PCA_D(1,Y==i),PCA_D(2,Y==i),PCA_D(3,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
    colormap(colorcube);
    hold on
end
legend('Basal','HER2','LA','LB','Normal-like','Normal');
hold on;
for i=1:N
    for j=2:size(IDX,2)
        plot3([PCA_D(1,IDX(i,1)) PCA_D(1,IDX(i,j))],[PCA_D(2,IDX(i,1)) PCA_D(2,IDX(i,j))],[PCA_D(3,IDX(i,1)) PCA_D(3,IDX(i,j))],'--k')
    end
end


grid;
end
function show_stru(X,Y,name)
FaceColor =  distinguishable_colors(max(Y));
flag=1;
switch name
    case 'PCA'
        [mapped_data,~]=compute_mapping(X','PCA',3);
        mapped_data=mapped_data';
    case 'LDA'
        [mapped_data,~]=compute_mapping([Y X'],'LDA',3);
        mapped_data=mapped_data';
    case 'DDR'
        params.maxIter = 20;
        params.eps = 1e-3;
        params.dim = 100;%rdim;
        params.lambda = 1 * size(X,2);% tree strength
        params.sigma = 0.1;%kernel width
        params.gamma = 4;% MSE
        [~, Z, ~, ~, ~] = DDRTree(X, params);
        mapped_data=Z;
        %         [mapped_data,~]=compute_mapping(Z','PCA',3);
        %         mapped_data=mapped_data';
    case 'PPT'
        params = [];
        params.maxIter = 100;
        params.lambda = 0.01;
        params.val=1;
        s = 0.01;
        params.bandwidth = 2 * s * s;
        % principal tree method
        [model,~] = principal_tree(X, params);
        % plot results
        stree = model.stree;
        mu = model.MU;
        Class = knnclassify(mu', X', Y,3);
        mapped_data=display_cancer([X mu],[Y;Class]);
        flag=0;
        if params.val==1
            PCA_D=mapped_data(:,length(Y)+1:end);
            hold on
            [m,n]=size(stree);
            for i=1:m
                for j=1:n
                    if stree(i,j)~=0
                        plot3([PCA_D(1,i) PCA_D(1,j)],[PCA_D(2,i) PCA_D(2,j)],[PCA_D(3,i) PCA_D(3,j)],'--k');
                    end
                end
            end
            hold off
        end
    case 'tSNE'
        [mapped_data,~]=compute_mapping(X','tSNE',3);
        mapped_data=mapped_data';
    case 'LLE'
        [mapped_data,~]=compute_mapping(X','tSNE',3);
        mapped_data=mapped_data';
end
if flag==1
    figure,
    for i=1:max(Y)
        plot3(mapped_data(1,Y==i),mapped_data(2,Y==i),mapped_data(3,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
        colormap(colorcube);
        hold on
    end
    legend('Basal','HER2','LA','LB','Normal-like','Normal');
    grid;
end
end
function mapped_data=display_cancer(X,Y,Name)
FaceColor =  distinguishable_colors(max(Y));
if size(X,1)>2
    type=3;
else
    type=2;
end
switch type
    case 3
        [mapped_data,~]=compute_mapping(X','PCA',3);
        mapped_data=mapped_data';
        figure,
        for i=1:max(Y)
            plot3(mapped_data(1,Y==i),mapped_data(2,Y==i),mapped_data(3,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
            colormap(colorcube);
            hold on
        end
    case 2
        [mapped_data,~]=compute_mapping(X','PCA',2);
        mapped_data=mapped_data';
        figure,
        for i=1:max(Y)
            plot(mapped_data(1,Y==i),mapped_data(2,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
            colormap(colorcube);
            hold on
        end
end
grid
if nargin<3
    legend('Basal','HER2','LA','LB','Normal-like','Normal');
else
    switch Name
        case 0
        case 1
            legend('Basal');
        case 2
            legend('Basal','Normal');
        case 4
            legend('1','2','3','4','Normal');
        case 5
            legend('Basal','HER2','LA','LB','Normal');
        case 6
            legend('Normal','6','7','8','9','10');
        case 9
            legend('1','2','3','4','HER2','LA','LB','Normal-like','Normal');
        otherwise
            legend('1','2','3','4','5','6','7','8','9','10');
    end
end
end
function heatmap_analy(X)
heatmap(X)
end
function cluster_analy(X,Y,type)
cg_s = clustergram(X, 'ColumnLabels', Y,...
    'Cluster',type,'Colormap',redbluecmap);
end
function plotSpanningTree(X,Y,dis)
if nargin<3
    dis='cityblock';
end
N=size(X,2);
Xw=(X)';
Fdis=pdist2(Xw,Xw,dis);
if ~strcmpi(dis,'cityblock')
    Fdis=Fdis.^2;
end
G=Fdis-diag(diag(Fdis));
[Tree, ~] = graphminspantree(sparse(G),'METHOD','Kruskal');
MST=full(Tree);
B=MST>0;
degree=sum(B,2);
display(['Max degree of the spanning tree ' num2str(max(degree))]);
display(['Mean degree of the spanning tree ' num2str(mean(degree))]);
[PCA_D,~]=compute_mapping(X','PCA',3);
PCA_D=PCA_D';
FaceColor =  distinguishable_colors(max(Y));
figure,
for i=1:max(Y)
    plot3(PCA_D(1,Y==i),PCA_D(2,Y==i),PCA_D(3,Y==i),'o','MarkerFaceColor',FaceColor(i,:));
    colormap(colorcube);
    hold on
end
legend('Basal','HER2','LA','LB','Normal-like','Normal');
hold on;
m=size(B,1);
for i=1:m
    for j=1:m
        if B(i,j)~=0
            plot3([PCA_D(1,i) PCA_D(1,j)],[PCA_D(2,i) PCA_D(2,j)],[PCA_D(3,i) PCA_D(3,j)],'--k');
        end
    end
end
grid;
hold off
end
function [W, Z, stree, Y, history] = DDRTree(X, params)
% X : DxN data matrix
% params.
%       maxIter : maximum iterations
%       eps     : relative objective difference
%       dim     : reduced dimension
%       lambda  : regularization parameter for inverse graph embedding
%       sigma   : bandwidth parameter
%       gamma   : regularization parameter for k-means
[~, N] = size(X);
% initialization
W = pca_projection(X * X', params.dim);
Z = W' * X;
if ~isfield(params,'ncenter')
    K = N;
    Y = Z(:,1:K);
else
    K = params.ncenter;
    [~, Y] = kmeans(Z',K);
    Y = Y';
end
% main loop
objs = [];
for iter=1:params.maxIter
    % Kruskal method to find optimal B
    distsqMU = sqdist(Y,Y);
    stree = graphminspantree(sparse(tril(distsqMU)),'Method','Kruskal');
    stree = stree + stree';
    B = stree ~= 0;
    L = diag( sum(B,2) ) - B;
    % compute R using mean-shift update rule
    distZY = sqdist(Z,Y);
    min_dist = repmat(min(distZY,[],2),1,K);
    tmp_distZY = distZY - min_dist;
    tmp_R = exp(-tmp_distZY ./ params.sigma);
    R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    Gamma = diag( sum(R) );
    % termination condition
    obj1 = - params.sigma * sum( log( sum( exp(- tmp_distZY./params.sigma) ,2) ) ...
        - min_dist(:,1) ./ params.sigma );
    objs(iter) = (norm(X-W*Z))^2 + params.lambda .* trace( Y * L * Y' )...
        + params.gamma * obj1;
    fprintf('iter=%d obj = %f\n',iter,objs(iter));
    history.W{iter} = W;
    history.Z{iter} = Z;
    history.Y{iter} = Y;
    history.stree{iter} = stree;
    history.R{iter} = R;
    if iter >1
        if abs(objs(iter) - objs(iter-1))/abs(objs(iter-1)) < params.eps
            break;
        end
    end
    % compute low dimension projection matrix
    tmp = R / ( ((params.gamma +1) / params.gamma) .* ( (params.lambda / params.gamma)  .* L + Gamma) - R' * R);
    Q = 1/(params.gamma +1) .* ( eye(N,N) + tmp * R' );
    C = X * Q;
    tmp1 = C * X';
    W = pca_projection( (tmp1 + tmp1')./2, params.dim);
    Z = W' * C;
    Y = Z * R / (params.lambda / params.gamma .* L + Gamma);
end
history.objs = objs;
end