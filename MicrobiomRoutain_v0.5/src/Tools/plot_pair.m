function plot_pair(D,IDX,index,MARK,dim)
if nargin<3
    index=(1:size(D,1))';
end
if nargin<4
    MARK='-or';
end
if nargin<5
    dim=3;
end
N=size(index,1);
if dim==2
    for i=1:N
        for k=1:size(IDX,2)
            plot([D(index(i),1) D(IDX(index(i),k),1)],[D(index(i),2) D(IDX(index(i),k),2)],MARK);
            hold on
        end
    end
    
else
    for i=1:N
        for k=1:size(IDX,2)
            plot3([D(index(i),1) D(IDX(index(i),k),1)],[D(index(i),2) D(IDX(index(i),k),2)],[D(index(i),3) D(IDX(index(i),k),3)],MARK);
            hold on
        end
    end

end
end