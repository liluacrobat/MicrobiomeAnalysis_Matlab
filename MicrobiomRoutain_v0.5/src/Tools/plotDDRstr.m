function [X,W,Z, stree,mu,tree_length]=plotDDRstr(X,Y,params)
if nargin<3
    params.maxIter = 20;
    params.eps = 1e-3;
    params.dim = size(X,1)/2;%rdim;
    params.lambda = 5 * size(X,2);% tree strength
    params.sigma = 0.001;%kernel width
    params.gamma = 4;% MSE
end
X = X - repmat( mean(X,2), [1 size(X,2)] );
[W, Z, stree,mu, ~,tree_length] = DDRTree(X, params);
mapped_data=Z;
if params.plot 
plotDDR(Z, mu, stree, Y);
end
end
function [W, Z, stree, Y, history,tree_length] = DDRTree(X, params)
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
    tree_length(iter)=trace( Y * L * Y' );
end
history.objs = objs;
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