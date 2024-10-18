function [alpha,order,tax_lasso,rpath] = CalLasso(x,y,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: m x n matriex with m features n samples
% y: 1 x n of sample labels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lasso_type (1:regression,2:classification)
if isfield(para,'lasso_type')
    lasso_type = para.lasso_type;
else
    lasso_type = 2;
end

if isfield(para,'lambda')
    lambda = para.lambda;
else
    lambda = 10;
end
if isfield(para,'tax')
    tax = para.tax;
else
    dim = size(x,1);
    tax = cell(d,1);
    for i=1:dim
        tax{i} = strcat('Feature',num2str(i));
    end
end
y_legend = para.y_legend;
if ~iscell(y_legend)
    tmp = cell(size(y_legend));
    for i=1:length(tmp)
        tmp{i} = num2str(y_legend(i));
    end
    y_legend = tmp;
end

rpath = para.rpath;
rpath = strcat(rpath,'_CMP');
[y_unique,~,y_rank] = unique(y);
y_cmp = y_legend(unique(y_unique));

[alpha] = HSICLasso(x,y_rank(:)',lasso_type,lambda);
order = find(alpha);
disp(['# of features selected by lasso: ' num2str(length(order))]);

if isempty(order)
    tax_lasso = [];
else
    tax_lasso = tax(order);
end
save(rpath,'alpha','order','tax_lasso','rpath','para');

if isfield(para,'plot')==0
    para.plot = 1;
end

if para.plot~=0
    measure_beta = x(order,:);
    measure_beta_title = 'euclidean';
    fileName = rpath;
    if size(x,1)>=2
        betaAnalysis(fileName,measure_beta,measure_beta_title,y_rank,y_cmp,2,1);
    end
end
end
