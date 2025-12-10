function [Error,CurveLength,mu,model,IDX]=plotPPTstr(X,Y,params)
IDX=[];
if nargin<3
    params = [];
    params.maxIter = 100;
    params.lambda = 0.01;
    params.val=1;
    s = 0.01;
    params.bandwidth = 2 * s * s;
    params.knn=1;
end
params.maxIter = 100;
% params.knn=0;
params.var=1;
knn=params.knn;
% principal tree method
switch params.var
    % %     case 1
    % [model,history] = principal_tree_var(X, params);
    %     case 2
    %         [model,history] = robust_principal_tree(X, params);
    %     case 3
    %         [model,history] = robust_principal_tree_noise(X, params);
    %         case 4
    %         [model,history] = robust_principal_tree_noise_p(X, params);
    %         case 5
    %         [model,history] = principal_graph(X, params);
    otherwise
        [model,history] = principal_tree(X, params);
end
% [model,history] = principalTreeMain(X, params);
% plot results
stree = model.stree;
mu = model.MU;
% Class = knnclassify(mu', X', Y,3);

if ~isfield(params,'val')
    params.val=1;
end
% knn=1;
mapped_data=plotPCA3D([X mu],[Y;ones(size(mu,2),1)*max(Y+1)],params.val);


if params.val==1
    PCA_X=mapped_data(:,1:length(Y));
    PCA_D=mapped_data(:,length(Y)+1:end);
    
    hold on
    [m,n]=size(stree);
    for i=1:m
        for j=1:n
            if stree(i,j)~=0
                plot3([PCA_D(1,i) PCA_D(1,j)],[PCA_D(2,i) PCA_D(2,j)],[PCA_D(3,i) PCA_D(3,j)],'-k','linewidth',3);
            end
        end
    end
    
    if knn>0
        IDX = knnsearch(mu',X','K',knn);
        for j=1:knn
            for i=1:size(X,2);
                plot3([PCA_X(1,i) PCA_D(1,IDX(i,j))],[PCA_X(2,i) PCA_D(2,IDX(i,j))],[PCA_X(3,i) PCA_D(3,IDX(i,j))],'--k','linewidth',1);
            end
        end
    end
    hold off
end
Error=history.mse(end);
CurveLength= history.length(end);
end
function mapped_data=plotPCA3D(X,Y,val)
U=sort(unique(Y));
str=cell(1,length(U));
FaceColor =  distinguishable_colors(length(U));
n=size(X,1);
if n<3
    X=[X;zeros(3-n,size(X,2)) ];
end

[mapped_data,~,power]=compute_mapping(X','PCA',3);
mapped_data=mapped_data';
if val==1
    figure,
    for i=1:length(U)
        plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',8,'MarkerEdgeColor','none');
        colormap(colorcube);
        hold on
        str{i}=num2str(U(i));
    end
    xlabel(['PC1 (' num2str(round(power(1)*1000)/10) '%)']);
    ylabel(['PC2 (' num2str(round(power(2)*1000)/10) '%)']);
    zlabel(['PC3 (' num2str(round(power(3)*1000)/10) '%)']);
    legend(str)
    grid
    boldify2
end
end
function [model,history] = principal_tree(X,params)

[D, N] = size(X);

% initialize MU, lmabda,
if ~isfield(params,'MU')
    K = N;
    MU = X;
else
    MU = params.MU;
    K = size(MU,2);
end

lambda = params.lambda * N; % scale the parameter by N

if ~isfield(params,'bandwidth')
    distsqX = sqdist(X,X);
    sigma = 0.01 * sum(sum(distsqX))/ (N*N);
else
    sigma = params.bandwidth;
end

for iter=1:params.maxIter
    
    % Kruskal method to find a spanning tree
    distsqMU = sqdist(MU,MU);
    stree = graphminspantree(sparse(tril(distsqMU)),'Method','Kruskal');
    stree = stree + stree';
    e = stree ~= 0;
    
    % NOTE: the raw implementation will have dividing zero problem
    %     % update data assignment matrix
    %     distMUX = sqdist(MU,X);
    %     tmp_R = exp(- distMUX'./sigma);
    %     R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    %
    %     % compute termination condition
    %     obj1 = - sigma * sum(log( sum(tmp_R',1) ) ) + 0.5 * lambda * sum(sum(full(stree)));
    %     % %-NOTE: probably log zero
    %     %  obj1 = sum(sum( R .* (distMUX' + sigma .* log (R)) + 0.5 * lambda .* full(stree)));
    %     history.objs(iter) = obj;
    %     % projected mean square error
    % %     projd = min(distMUX,[],1);
    % %     mse = mean(projd);
    % %     history.mse(iter) = mse;
    %     mse1 = sum(sum( R .* distMUX'));
    %     if isnan(mse)
    %         fprintf('for check!\n');
    %     end
    %     history.mse(iter) = mse;
    
    
    history.mu{iter} = MU;
    history.stree{iter} = stree;
    
    % update data assignment matrix
    distMUX = sqdist(MU,X);
    min_dist = repmat(min(distMUX,[],1),K,1);
    tmp_distMUX = distMUX - min_dist;
    tmp_R = exp(-tmp_distMUX' ./ sigma);
    R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    
    % compute objective function
    obj1 = - sigma * sum( log( sum( exp(- tmp_distMUX./sigma) ,1) ) ...
        - min_dist(1,:) ./ sigma );
    reg = sum(sum(full(stree)));
    obj = (obj1 + 0.5 * lambda * reg)/N;
    history.objs(iter) = obj;
    
    % projected mean square error
    %     mse = sum(sum( R .* distMUX'));
    %     mse = obj1 + sigma *N* log(K);
    
    projd = min(distMUX,[],1);
    mse = mean(projd);
    
    history.mse(iter) = mse;
    
    % length of the structure
    history.length(iter)=reg;
    
    % terminate condition
    fprintf('iter=%d, obj=%f, mse=%f, len=%f\n',iter,obj,mse,reg);
    if iter > 1
        if abs((obj - old_obj)/old_obj) < 1e-3
            break;
        end
    end
    
    % compute the means
    L = diag(sum(e,1)) - e;
    MU = X * R / ( lambda * L + diag(sum(R,1)) );
    
    old_obj = obj;
end

model.MU = MU;
model.stree = stree;
model.sigma = sigma;
model.lambda = params.lambda;
end
function iPT_debug(Xw,Fw)
plotPCA3D([Xw Fw],[ones(size(Xw,2),1);ones(size(Fw,2),1)*2],1);
legend('Sample','Tree');
end