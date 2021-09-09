function [tstat,para]=StaTest(X,Y,method)
tstat=zeros(size(X,1),size(Y,2));
para=[];
switch method
    case 'kruwallis'
        for i=1:size(X,1)
            [p,tbl,stats] = kruskalwallis(X(i,:),Y,'off');
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
        [rho,tstat]=corr(X',Y(:),'type','Spearman');
        para=rho;
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