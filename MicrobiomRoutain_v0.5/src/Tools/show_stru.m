function Fin=show_stru(X,Y,name,Para)
Fin=[];
FaceColor =  distinguishable_colors(max(Y));
flag=1;
switch name
    case 'KNN'
        k=Para;
        [IDX]=plot_knn(X,Y,k);
    case 'spanning tree'
        dis=Para;
        plotSpanningTree(X,Y,dis);
    case 'heatmap'
        heatmap_analy(X,Para);
    case 'cluster'
        type=Para;
        Fin=cluster_analy(X,Y,type);
    otherwise
        switch name
            case 'PCA'
                [mapped_data,~]=compute_mapping(X','PCA',4);
                mapped_data=mapped_data';
            case 'LDA'
                [mapped_data,~]=compute_mapping([Y X'],'LDA',3);
                mapped_data=mapped_data';
            case 'DDR'
                if nargin<4
                params.maxIter = 20;
                params.eps = 1e-3;
                params.dim = size(X,1)/2;%rdim;
                params.lambda = 5 * size(X,2);% tree strength
                params.sigma = 0.001;%kernel width
                params.gamma = 4;% MSE
                else
                    params=Para;
                end
                [~, Z, stree,mu, ~] = DDRTree(X, params);
                mapped_data=Z;
                plotDDR(Z, mu, stree, Y);
                flag=0;
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
            case 'isomap'
                [mapped_data,mapping]=compute_mapping(X','Isomap',3,Para);
                mapped_data=mapped_data';
                flag=0;
                figure,
                plot3(mapped_data(1,:),mapped_data(2,:),mapped_data(3,:),'o','MarkerFaceColor',FaceColor(1,:));
                colormap(colorcube);
                
        end
        if flag==1
            figure,
            U=unique(Y);
            str=cell(1,length(U));
            for i=1:length(U)
                plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:));
                colormap(colorcube);
                hold on
                str{i}=num2str(U(i));
            end
            legend(str);
            grid;
        end
end
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
function heatmap_analy(X,Para)
xl=Para.x;
yl=Para.y;
heatmap(X,xl,yl);
end
function cg_s=cluster_analy(X,Y,type)
% cg_s = clustergram(X, 'ColumnLabels', Y,...
%     'Cluster',type,'Colormap',redbluecmap,'Symmetric',false);
cg_s = clustergram(X, ...
    'Cluster',type,'Colormap',redbluecmap,'Symmetric',false);
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
function plotDDR(X, mu, stree, PAM50)
% draw cancer data
% line_config ={'>','s','o','d','^','p'};
% names = {'1','2','3','4','5','6'};
% names = {'Basal','Her2+','luminal A', 'luminal B', 'normal-like','normal'};
mpdc6 = distinguishable_colors(length(unique(PAM50)));
U=unique(PAM50);
U=sort(U);
numL = length(U);

figure;
for i=1:numL
    idx = find(PAM50==U(i));
    names{i}=num2str(U(i));
    hp(i)=plot3(X(1,idx), X(2,idx), X(3,idx),'o',...
        'MarkerFaceColor',mpdc6(i,:));
    hold on;
end
grid on;

for n = 1:size(stree)
    index = find(stree(n,:)>0);
    for m = 1:length(index)
        hline=plot3([mu(1,n), mu(1,index(m))],[mu(2,n), mu(2,index(m))],...
            [mu(3,n), mu(3,index(m))],'-k','LineWidth',3);
    end
end

% set(gca, 'FontSize',16);
legend([hp hline],[ names(1:numL), 'tree structure']);
end