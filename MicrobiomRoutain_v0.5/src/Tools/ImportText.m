function R=ImportText(name,start,style,skip)
id=fopen(name,'r');
if nargin<2
    start=1;
end
if nargin<3
    style='num';
end
if start>1
    for i=1:start-1
        line=fgetl(id);
    end
end
if nargin<4
    skip=0;
end
n=0;
while ~feof(id)
    n=n+1;
    line=fgetl(id);
    S=strsplit(line,'\t');
    if  n==1
        if skip==0
            R.colname=S(2:end);
            continue
        else
            n=n+1;
        end
    end
    R.rowname{n-1,1}=S{1};
    num=length(S)-1;
    if strcmpi(style,'txt')
        for i=1:num
            R.Data{n-1,i}=S{i+1};
        end
    else
        
        for i=1:num
            R.Data(n-1,i)=str2num(S{i+1});
        end
    end
    
end
end