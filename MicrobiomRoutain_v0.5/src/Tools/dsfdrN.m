function [reject,pvals,result]=dsfdrN(data,labels,method,alpha,fdr_method,numperm)
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
rng(100)
if nargin<6
    numperm=1000;
end
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

if isempty(fdr_method)
    fdr_method='dsfdr';
end

if isempty(method)
    method='meandiff';
end
[numbact,S]=size(data);

labels=labels(:);
switch method
    case 'kruwallis'
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            display(['Perm: ' num2str(i)]);
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
    case 'ranksum'
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            display(['Perm: ' num2str(i)]);
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
    case 'ttest'
        tstat=StaTest(data,labels,method);
        t=abs(tstat);
        u=zeros(numbact,numperm);
        for i=1:numperm
            display(['Perm: ' num2str(i)]);
            rlabels=labels(randperm(S));
            rt=StaTest(data,rlabels,method);
            u(:,i)=rt;
        end
        
    case 'spearman'
        [tstat,rho]=StaTest(data,labels,method);
        %         data=para_spearman.x;
        %         labels=para_spearman.y;
        t = abs(tstat);
        permlabels = zeros([length(labels), numperm]);
        for cperm =1: numperm
            display(['Perm: ' num2str(cperm)]);
            rnd_num = randperm(length(labels));
            rlabels = labels(rnd_num);
            permlabels(:, cperm) = rlabels(:);
        end
        %         u = abs(data, permlabels);
        u=StaTest(data,permlabels,method);
        result.rho=zeros(size(idx_s))*nan;
        result.rho(idx_s==0)=rho;
end

pvals=tstat;
pvals_u=u;



% calculate FDR
% pvals=round(pvals/10^(-8))*10^(-8);
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
    if fdr <alpha&&foundit ==0
        realcp = cp;
        foundit = 1;
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
result.tstat=zeros(size(idx_s))*nan;
result.tstat(idx_s==0)=tstat;
result.s=idx_s;

result.pvals=zeros(size(idx_s))*nan;
result.pvals(idx_s==0)=pvals;

result.pvals_u=zeros(size(idx_s,1),numperm)*nan;
result.pvals_u(idx_s==0,:)=pvals_u;

result.realcp=realcp;
result.allfdr=allfdr;
result.allt=allt;
result.Q=zeros(size(idx_s))*nan;
result.Q(idx_s==0)=adjp;
end

function [tstat,para]=StaTest(X,Y,method)
tstat=zeros(size(X,1),size(Y,2));
para=[];
U=unique(Y);
switch method
    case 'kruwallis'
        for i=1:size(X,1)
            [p,tbl,stats] = kruskalwallis(X(i,:),Y,'off');
            %             p=ranksum(X(i,Y==U(1)),X(i,Y==U(2)));
            tstat(i)=p;%tbl{2,5};
        end
    case 'ranksum'
        for i=1:size(X,1)
            p=ranksum(X(i,Y==U(1)),X(i,Y==U(2)));
            tstat(i)=p;%tbl{2,5};
        end
    case 'ttest'
        U=unique(Y);
        for i=1:size(X,1)
            [~,p] = ttest2(X(i,Y==U(1)),X(i,Y==U(2)));
            tstat(i)=p;%tbl{2,5};
        end
    case 'spearman'
        %         Y=transform_rankdata(Y);
        %         Y=Y-mean(Y);
        %         X=transform_rankdata(X);
        %         meanval=mean(X,2);
        %         X=X-repmat(meanval,1,size(X,2));
        %         tstat=X*Y(:);
        %         para.x=X;
        %         para.x=Y;
        [para,tstat]=corr(X',Y,'type','Spearman');
        
    case 'nonzerospearman'
        for i=1:size(X,1)
            idx0=find(X(i,:)~=0);
            Xt=X(i,idx0);
            Yt=Yt(idx0);
            
            Xt=transform_rankdata(Xt);
            Yt=transform_rankdata(Yt);
            Yt=Yt-mean(Yt);
            
            meanval=mean(Xt);
            Xt=Xt-meanval;
            tstat(i)=Xt*Yt(:);
        end
end
end

