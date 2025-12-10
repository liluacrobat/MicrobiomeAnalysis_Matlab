function [model,history] = principal_tree_var(X,params)

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

old_length = 1;
for iter=1:params.maxIter
    iter
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
    distMUX=distMUX';
    min_dist = repmat(min(distMUX,[],1),K,1);
    tmp_distMUX = distMUX - min_dist;
    tmp_R = exp(-tmp_distMUX' ./ sigma);
    R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    distMUX=distMUX';
    R1=R';
    min_dist = repmat(min(distMUX,[],1),K,1);
    tmp_distMUX = distMUX - min_dist;
    tmp_R = exp(-tmp_distMUX' ./ sigma);
    R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    R2=R;
    R=(params.bal*R1+(1-params.bal)*R2);
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
    length = history.length(end);
    
    % terminate condition
    %fprintf('iter=%d, obj=%f, mse=%f, len=%f\n',iter,obj,mse,reg);
    if iter > 1
        %         if abs((obj - old_obj)/old_obj) < 1e-10
        %             break;
        %         end
        if abs((length - old_length)/old_length) < 1e-20
            break;
        end
    end
    
    % compute the means
    L = diag(sum(e,1)) - e;
    MU = X * R / ( lambda * L + diag(sum(R,1)) ); 
    
    old_obj = obj;
    old_length = length;
end

model.MU = MU;
model.stree = stree;
model.sigma = sigma;
model.lambda = params.lambda;

%% ========= The end of the code ==================


