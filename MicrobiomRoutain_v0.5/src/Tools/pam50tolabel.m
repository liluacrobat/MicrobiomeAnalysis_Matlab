function pam50tolabel(pam50)
n=length(pam50);
Y=zeros(n,1);
for i=1:n
    str=pam50{i};
    switch str
        case 'Basal'
            Y(i)=1;
        case 'Her2'
            Y(i)=2;
        case 'LumA'
            Y(i)=3;
        case 'LumB'
            Y(i)=4;
        case 'Normal'
            Y(i)=5;
        case 'NT'
            Y(i)=0;
        otherwise
            Y(i)=nan;
    end
end
end