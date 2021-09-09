function [reject,pvals,result]=dsfdr(data,labels,method,alpha,fdr_method)
%%
% calculate the Discrete FDR for the data
%
%     input:
%     data : N x S numpy array
%         each column is a sample (S total), each row an OTU (N total)
%     labels : a 1d numpy array (length S)
%         the labels of each sample (same order as data) with the group
%         (0/1 if binary, 0-G-1 if G groups, or numeric values for correlation)
%
%
%     transform_type : str or None
%         transformation to apply to the data before caluculating
%         the test statistic
%         'rankdata' : rank transfrom each OTU reads
%         'log2data' : calculate log2 for each OTU using minimal cutoff of 2
%         'normdata' : normalize the data to constant sum per samples
%         'binarydata' : convert to binary absence/presence
%         'clrdata' : clr transformation of data (after replacing 0 with 1)
%          None : no transformation to perform
%
%     method : str or function
%         the method to use for calculating test statistics:
%         'meandiff' : mean(A)-mean(B) (binary)
%         'mannwhitney' : mann-whitney u-test (binary)
%         'kruwallis' : kruskal-wallis test (multiple groups)
%         'stdmeandiff' : (mean(A)-mean(B))/(std(A)+std(B)) (binary)
%         'spearman' : spearman correlation (numeric)
%         'pearson' : pearson correlation (numeric)
%         'nonzerospearman' : spearman correlation only non-zero entries
%                             (numeric)
%         'nonzeropearson' : pearson correlation only non-zero entries (numeric)
%         function : use this function to calculate the test statistic
%         (input is data,labels, output is array of float)
%
%     alpha : float
%         the desired FDR control level
%     numperm : int
%         number of permutations to perform
%
%     fdr_method : str
%         the FDR procedure to determine significant bacteria
%         'dsfdr' : discrete FDR method
%         'bhfdr' : Benjamini-Hochberg FDR method
%         'byfdr' : Benjamini-Yekutielli FDR method
%         'filterBH' : Benjamini-Hochberg FDR method with filtering
%         'gilbertBH' : Benjamini-Hochberg FDR method with Gilbert (2005) pre-filtering
%
%     output:
%     reject : np array of bool (length N)
%         True for OTUs where the null hypothesis is rejected
%     tstat : np array of float (length N)
%         the test statistic value for each OTU (for effect size)
%     pvals : np array of float (length N)
%         the p-value for each OTU
%%
rng default
idx_s=zeros(size(data,1),1);
for i=1:size(data,1)
    if length(unique(data(i,:)))==1
        idx_s(i)=1;
    end
end
data=data(idx_s==0,:);
orig_numbact=size(data,1);
filtered_order =1:orig_numbact;
if strcmp(fdr_method,'filterBH')==1
    index=[];
    U=unique(labels);
    temp=labels;
    labels(temp==U(1))=0;
    labels(temp==U(2))=1;
    n0=sum(labels==0);
    n1=sum(labels==1);
    for i=1:orig_numbact
        nonzeros = length(find(data(i, :)~=0));
        if nonzeros<min(n0,n1)
            pval_min =(length(combnk(n0, nonzeros))+length(combnk(n1, nonzeros)))...
                /length(combnk(n0+n1, nonzeros));
            if pval_min<=alpha
                index=[index i];
            end
        else
            index=[index i];
        end
    end
    data=data(index,:);
    filtered_order=filtered_order(index);
else
    if strcmp(fdr_method,'gilbertBH')==1
        % caluclate the Gilbert alpha* per feature (minimal ibtainable p-value)
        alpha_star = [];
        U=unique(labels);
        temp=labels;
        labels(temp==U(1))=0;
        labels(temp==U(2))=1;
        n0=sum(labels==0);
        n1=sum(labels==1);
        for i=1:orig_numbact
            % test if all values are identical, max p-val is 1 (need to filter)
            if length(unique(data(i,:)))==1
                alpha_star=[alpha_star 1];
                continue
            end
            cdat = sort(data(i,:)); % sort in acending order
            rdat = sort(data(i,:),'descend'); % sort in decending order
            
            [p1,tbl1,stats1] = kruskalwallis(cdat(:)',[zeros(1,n0) ones(1,n1)],'off');
            [p2,tbl2,stats2] = kruskalwallis(rdat(:)',[zeros(1,n0) ones(1,n1)],'off');
            s1=tbl1{2,5};
            s2=tbl2{2,5};
            alpha_star=[alpha_star min(p1,p2)];
            %find the smallest K which is big enough for Bonferoni (that's how it's done in Gilbert)
            for ck = 1:orig_numbact
                num_ok = sum(alpha_star < alpha / ck);
                if num_ok <= ck
                    break
                end
                
                
            end
            %and keep only the features which match it
            index = (alpha_star < alpha / ck);
            data = data(index, :);
            filtered_order = filtered_order(index);
        end
        
    end
end

% transform the data

if isempty(alpha)
    alpha=0.05;
end
numperm=1000;
if isempty(fdr_method)
    fdr_method='dsfdr';
end

if isempty(method)
    method='meandiff';
end
[numbact,S]=size(data);


switch method
    
    case 'spearman'
        [tstat,para_spearman]=StaTest(data,labels,method);
        data=para_spearman.x;
        labels=para_spearman.y;
        t = abs(tstat);
        permlabels = zeros([length(labels), numperm]);
        for cperm =1: numperm
            rlabels = randperm(labels);
            permlabels(:, cperm) = rlabels(:);
        end
        u = abs(data*permlabels);
    otherwise
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
        for i=1:numbact
            R=isclose(t(i),u(i,:));
            u(i,R==1)=t(i);
        end
end

% calculate permutation p-vals
pvals = zeros(numbact,1);
pvals_u = zeros([numbact, numperm]);
% pseudo p-values for permutated test statistic u
for i=1:numbact
    allstat = [t(i), u(i,:)];
    stat_rank=Rankings(allstat,'competition1');
    allstat = 1 - ((stat_rank - 1) / length(allstat));
    % assign ranks to t from biggest as 1
    pvals(i) = allstat(1);
    pvals_u(i,:) = allstat(2:end);
end
% calculate FDR
pvals_unique = unique(pvals);
sortp = sort(pvals_unique,'descend');
% find a data-dependent threshold for the p-value
foundit = 0;
allfdr = [];
allt = [];
realcp=2;
for i=1:length(sortp)
    cp=sortp(i);
    realnum = sum(pvals <= cp);
    fdr = (realnum + length(find(pvals_u <= cp))) / ...
        (realnum * (numperm + 1));
    allfdr=[allfdr fdr];
    allt=[allt cp];
    if fdr <alpha
        realcp = cp;
        foundit = 1;
        break
    end
end

if foundit ==0
    reject=zeros(numbact,1);
else
    reject=pvals <= realcp;
end
adjp=pvals;
for i=1:length(allt)
    adjp(pvals==allt(i))=allfdr(i);
end
result.s=idx_s;
result.tstat=tstat;
result.pvals_u=pvals_u;
result.realcp=realcp;
result.allfdr=allfdr;
result.allt=allt;
result.Q=adjp;
end



