function fun_genNetwork(fname,r,v1,v2,thr,symt)
if nargin<5
    symt=1;
end
if nargin<4
    thr=0;
end
ff1 = fopen(strcat(fname,'_network.txt'),'w');
fprintf(ff1,'source\ttarget\tinteraction\tdirected\tweight\tw_sign\tedge_type\n');
[m,n]=size(r);
sn=0;
if symt==1
    for i=1:m
        for j=i+1:n
            if abs(r(i,j))>thr
                sn=sn+1;
                fprintf(ff1,'%s\t%s\tpp\tfalse\t%f\t%d\tg2c\n',v1{i},v2{j},r(i,j),sign(r(i,j)));
            end
        end
    end
else
    for i=1:m
        for j=1:n
            if abs(r(i,j))>thr
                sn=sn+1;
                fprintf(ff1,'%s\t%s\tpp\tfalse\t%f\t%d\tg2c\n',v1{i},v2{j},r(i,j),sign(r(i,j)));
            end
        end
    end

end

ff2 = fopen(strcat(fname,'_network_node.txt'),'w');
fprintf(ff2,'shared name\tdep\tnode_type\n');
for i=1:m
    fprintf(ff2,'%s\t%s\tMeta\n',v1{i},v1{i});
end
if symt~=1
    for i=1:n
        fprintf(ff2,'%s\t%s\tMeta\n',v2{i},v2{i});
    end
end
disp(sn);
end