function [ACC, w_history, Predict_value]=CrossVail_LOGO(X, label, Para)
%% ==========================================================================
% CrossVail_LOGO: Determined the parameters of feature selection using LOGO through
% crossvalidation
%
%--------------------------------------------------------------------------
% INPUT:
%     X:  training data: [x1,x2,...xn] Each column is an observation
%     label:  class label = {1,2,...,C}
%     Para:  parameters. More specifically,
%            Para.distance: distance (eg. block, euclidean)
%            Para.sigma: kernel width (k nearest neighbor)
%            Para.lambda: regulariztion parameter
%            Para.kernel: kernel
%            Para.maxit: maximum number of iteration
%            Para.folds: number of folds
%            Para.soft: 1: exponential kernel; 0: limited kernel
%            Para.plot: 1: plot of the learning process; 0: do not plot
%
% OUTPUT:
%     ACC: crossvalidation accuracy of each fold
%     w_history: feature weight of each fold
%     Predict_value: prediction value
%--------------------------------------------------------------------------
% by Lu Li
% update history: 7/23/2018
%% ==========================================================================
% partitioning samples into subsets
fold_id = selectFOLD(label, Para.folds);
nk = length(unique(label)); % number of label class
Para.plotfigure = 0;
%% initializations
ACC = zeros(1,max(fold_id));
w_history = zeros(size(X,1),10);
Predict_value = zeros(size(label));
%% apply crossvalidation
for k=1:max(fold_id)
    display(['Forld ' num2str(k)]);
    Para.plotfigure = 0;
    % partitioning samples into training and testing
    train_patterns = X(:,fold_id~=k);
    train_targets = label(fold_id~=k);
    test_patterns = X(:,fold_id==k);
    test_targets = label(fold_id==k);
    Para_train = Para;
    %% feature selection using LOGO
    [Weight,~,~] = Logo_kernel(train_patterns, train_targets, Para_train);
    Weight = Weight(:);
    w_history(:,k) = Weight;
    %% classification
    index_all = cell(1,nk);
    N = zeros(1,nk);
    patterns_nk = cell(1,nk);
    for i=1:nk
        index_all{i} = find(train_targets == i);
        N(i) = length(index_all{i});
        patterns_nk{i} = train_patterns(:, index_all{i});
    end
    Weight = Weight(:);
    Prho = zeros(size(test_patterns,2),nk);
    % predict based on the training data
    for n = 1:size(test_patterns,2)
        test = test_patterns(:,n);
        for i=1:nk
            % calculate the distance to samples with same label
            temp_H = abs(patterns_nk{i}-test*ones(1,N(i)));
            dist_H = (Weight)'*temp_H;
            % calculate the distance to samples with different label
            ID_M = find(train_targets~=i);
            ID_n = length(ID_M);
            patterns_M = train_patterns(:,ID_M);
            temp_M = abs(patterns_M-test*ones(1,ID_n));
            dist_M = (Weight)'*temp_M;
            
            prob_H = kernel_fun(dist_H,Para_train);
            prob_M = kernel_fun(dist_M,Para_train);
            
            rho = sum(dist_M.*prob_M)-sum(dist_H.*prob_H);
            Prho(n,i) = 1/(1+exp(-rho));
        end
    end
    % compute testing error
    Prediction = zeros(size(Prho,1),1);
    for i = 1:length(Prediction)
        [~,Prediction(i)] = max(Prho(i,:),[],2);
    end
    Predict_value(fold_id==k) = Prediction;
    Test_Error = length(find(Prediction(:)~=test_targets(:)))/length(test_targets);
    Test_Result = Test_Error;
    R = (Test_Result*100);
    ACC(k) = 100-R;
end
end

function fold_id = selectFOLD(label, folds)
%% ==========================================================================
%
% selectFOLD: partitioning samples into subsets
%
%--------------------------------------------------------------------------
%INPUT:
%     label:  class label = {1,2,...,C}
%     seed: random seed
%     folds: number of folds fro crossvalidation
%
%OUTPUT:
%     fold_id: partition of folds
%
%--------------------------------------------------------------------------
% by Lu Li
% update history: 7/23/2018
%% ==========================================================================
idx = 1:length(label);
rng(98,'twister');
nQ = length(unique(label));
l_G = zeros(1,nQ);
perm_G = cell(1,nQ);
step_G = zeros(1,nQ);
idx_G = cell(1,nQ);
% random select samples for each class
for i = 1:nQ
    l_G(i) = length(find(label==i));
    perm_G{i} = randperm(l_G(i));
    step_G(i) = round(l_G(i)/folds);
    idx_G{i} = idx(label==i);
end
% assign samples for each fold
fold_id = ones(length(label),1)*folds;
for i = 1:(folds-1)
    for j = 1:length(unique(label))
        idx_B = idx_G{j};
        perm_B = perm_G{j};
        l_B = l_G(j);
        fold_id(idx_B(perm_B(round(l_B/folds*(i-1))+1:round(l_B/folds*i)))) = i;
    end
end
%% ==================End of the code===================================
end
