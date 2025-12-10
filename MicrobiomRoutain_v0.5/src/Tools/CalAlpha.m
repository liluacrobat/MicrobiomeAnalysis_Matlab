function [rich,even,rich2]=CalAlpha(X,para)
rng(100)
[N,M]=size(X);
switch para.rare
    case 1
        num=para.num;
        richT=zeros(M,para.rep);
        evenT=zeros(M,para.rep);
        richT2=zeros(M,para.rep);
        C=cumsum(X,1);
        for i=1:para.rep
            disp(['iteration ' num2str(i)]);
            X_s=zeros(size(X));
            for j=1:size(X,2)
                if num<=C(end,j)
                    temp=X(:,j);
                    for k=1:num
                        nz=find(temp>0);
                        sel=randperm(length(nz),1);
                        temp(nz(sel))=temp(nz(sel))-1;
                        X_s(nz(sel),j)=X_s(nz(sel),j)+1;
                    end
                else
                    X_s(:,j)=0;
                end
            end
            richT(:,i)=ObserveOTU(X_s);
            richT2(:,i)=CalChao1(X_s);
            [H,~]=index_SaW(X_s);
            evenT(:,i)=H;
        end
        rich=mean(richT,2);
        even=mean(evenT,2);
        rich2=mean(richT2,2);
        rich(num>C(end,:))=nan;
        rich2(num>C(end,:))=nan;
        even(num>C(end,:))=nan;
    case 2
        num=para.num;
        richT=zeros(M,para.rep);
        richT2=zeros(M,para.rep);
        evenT=zeros(M,para.rep);
        C=cumsum(X,1);
        Xfull=cell(1,M);
        for j=1:size(X,2)
            temp=zeros(C(end,j),1);
            temp(1:C(1,j))=1;
            for k=2:N
                temp(C(k-1,j)+1:C(k,j))=k;
            end
            Xfull{j}=temp;
        end
        
        for i=1:para.rep
            disp(['iteration ' num2str(i)]);
            X_s=zeros(size(X));
            for j=1:size(X,2)
                tic
                if num<=C(end,j)
                    sel=randperm(C(end,j),num);
                    idx=zeros(C(end,j),1);
                    idx(sel)=1;
                    temp=Xfull{j};
                    temp(idx==0)=0;
                    u=setdiff(unique(temp),0);
                    for k=1:length(u)
                        X_s(k,j)=length(find(temp==u(k)));
                    end
                else
                    X_s(:,j)=0;
                end
                toc
            end
            richT(:,i)=ObserveOTU(X_s);
            richT2(:,i)=CalChao1(X_s);
            [H,~]=index_SaW(X_s);
            evenT(:,i)=H;
        end
        rich=mean(richT,2);
        rich2=mean(richT2,2);
        even=mean(evenT,2);
        rich(num>C(end,:))=nan;
        rich2(num>C(end,:))=nan;
        even(num>C(end,:))=nan;
    otherwise
        rich=ObserveOTU(X);
        rich2=CalChao1(X);
        [H,~]=index_SaW(X);
        even=H;
end
end
function A=ObserveOTU(D)
A=zeros(1,size(D,2));
for i=1:length(A)
    s=length(find(D(:,i)>0));
    A(i)=s;
end
end
