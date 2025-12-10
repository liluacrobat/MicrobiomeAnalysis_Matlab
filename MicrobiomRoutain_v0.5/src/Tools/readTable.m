function [ SAMPLE_ID,OTU_ID,OTU_TABLE,TAXONOMY]=readTable(file_name,flag)
f_id=fopen(file_name,'r');
n=0;
if nargin<2
flag=0;
end
if flag==0
while ~feof(f_id)
    n=n+1;
    display(n);
    line=fgetl(f_id);
    S=strsplit(line,'\t');
    if  n==1
        SAMPLE_ID=S(1,2:end);
        OTU_TABLE=zeros(1,length(SAMPLE_ID));
    else
        OTU_ID{n-1}=S{1};
        TAXONOMY{n-1}=S{1};
        OTU_TABLE(n-1,:)=cell2Tnum(S(2:end));
    end
end
else
while ~feof(f_id)
    n=n+1;
    display(n);
    line=fgetl(f_id);
    S=strsplit(line,'\t');
    if  n==1
        SAMPLE_ID=S(1,2:end-1);
        OTU_TABLE=zeros(1,length(SAMPLE_ID));
    else
        OTU_ID{n-1}=S{1};
        TAXONOMY{n-1}=S{end};
        OTU_TABLE(n-1,:)=cell2Tnum(S(2:end-1));
    end
end
end
end
function Y=cell2Tnum(X)
Y=zeros(size(X));
for i=1:size(Y,1)
    for j=1:size(Y,2)
        if strcmp(X{i,j},'not collected')==1 || strcmp(X{i,j},'not applicable')==1
            Y(i,j)=nan;
        else
            Y(i,j)=str2num(X{i,j});
        end
    end
end
end
