function [observed_otu, shannon_index, chao1_index,simposon,Result] = CalAlphaDiversity(...
    X,para)
%% calculate alpha diversity
rng default
[N,M]=size(X);
Result = [];
if isfield(para,'num')
    para.rare = 2;
end
switch para.rare
    case 1
        %         num = para.num;
        %         X = round(CalRel(X)*num);
        observed_otu = ObserveOTU(X);
        chao1_index = CalChao1(X);
        [H,~] = index_SaW(X);
        shannon_index = H;
        simposon = CalSimpson(X);
        Result.otu = observed_otu;
        Result.shannon = shannon_index;
        Result.chao1 = chao1_index;
        Result.simposon = simposon;
    case 2
        num = para.num;
        observed_otu_t = zeros(M,para.rep);
        chao1_index_t = zeros(M,para.rep);
        shannon_index_t = zeros(M,para.rep);
        simposon_index_t = zeros(M,para.rep);
        C = round(cumsum(X,1));
        Xfull = cell(1,M);
        for j = 1:size(X,2)
            temp = zeros(C(end,j),1);
            temp(1:C(1,j)) = 1;
            for k = 2:N
                temp(C(k-1,j)+1:C(k,j)) = k;
            end
            Xfull{j} = temp;
        end
        
        for i = 1:para.rep
            disp(['iteration ' num2str(i)]);
            X_s = zeros(size(X));
            for j = 1:size(X,2)
                %                 tic
                if num <= C(end,j)
                    sel = randperm(C(end,j),num);
                    idx = zeros(C(end,j),1);
                    idx(sel) = 1;
                    temp = Xfull{j};
                    temp(idx==0) = 0;
                    u = setdiff(unique(temp),0);
                    for k = 1:length(u)
                        X_s(u(k),j) = length(find(temp==u(k)));
                    end
                else
                    num1 = C(end,j);
                    sel = randperm(C(end,j),num1);
                    idx = zeros(C(end,j),1);
                    idx(sel) = 1;
                    temp = Xfull{j};
                    temp(idx==0) = 0;
                    u = setdiff(unique(temp),0);
                    for k = 1:length(u)
                        X_s(u(k),j) = length(find(temp==u(k)));
                    end
                end
                %                 toc
            end
            observed_otu_t(:,i) = ObserveOTU(X_s);
            chao1_index_t(:,i) = CalChao1(X_s);
            [H,~] = index_SaW(X_s);
            shannon_index_t(:,i) = H;
            simposon_index_t(:,i) = CalSimpson(X_s);
        end
        Result.otu = observed_otu_t;
        Result.shannon = shannon_index_t;
        Result.chao1 = chao1_index_t;
        Result.simposon = simposon_index_t;
        observed_otu = mean(observed_otu_t,2);
        chao1_index = mean(chao1_index_t,2);
        shannon_index = mean(shannon_index_t,2);
        simposon = mean(simposon_index_t,2);
        observed_otu(num>C(end,:)) = nan;
        chao1_index(num>C(end,:)) = nan;
        shannon_index(num>C(end,:)) = nan;
        simposon(num>C(end,:)) = nan;
    otherwise
        observed_otu = ObserveOTU(X);
        chao1_index = CalChao1(X);
        [H,~] = index_SaW(X);
        shannon_index = H;
        simposon = CalSimpson(X);
end
end
function A = ObserveOTU(D)
A = zeros(1,size(D,2));
for i = 1:length(A)
    s = length(find(D(:,i)>0));
    A(i) = s;
end
end
function rich=CalChao1(X)
% PATH = getenv('PATH');
% setenv('PATH', [PATH ':/Library/Frameworks/R.framework/Resources/bin']);
save('Chao1_temp.mat','X');
system('R CMD BATCH EsChao1.R');
load('Chao1_index','rich');
delete('Chao1_temp.mat');
delete('Chao1_index.mat');
end
function rich=CalSimpson(X)
% PATH = getenv('PATH');
% setenv('PATH', [PATH ':/Library/Frameworks/R.framework/Resources/bin']);
save('Simpson_temp.mat','X');
system('R CMD BATCH EsSimpson.R');
load('Simpson_index','rich');
delete('Simpson_temp.mat');
delete('Simpson_index.mat');
end


