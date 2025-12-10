function genNet4Cyto(edges,node,para)
[m,n] = size(edges);
if isfield(para,'node2')
    node2 = para.node2;
else
    node2 = node;
end
if isfield(para,'thr')
    thr = para.thr;
else
    thr = 0.3;
end
if isfield(para,'dir')
    f_direct = para.dir;
else
    f_direct = 'false';
end
fid = fopen(strcat(para.file,'.txt'),'w');
fprintf(fid,'source\ttarget\tinteraction\tdirected\tweight\tedgedir\teffectsize\tsource_node\tcluster\n');
for i=1:m
    for j=1:n
        if abs(edges(i,j))>=thr
            fprintf(fid,'%s\t%s\tpp\t%s\t%f\t%d\t%f\tYes\t%f\n',...
                node{i},node2{j},f_direct,edges(i,j),sign(edges(i,j)),abs(edges(i,j)),para.T(j));
        end
    end
end

ee = abs(edges)>=thr;

if isfield(para,'node2')
    n1 = sum(sum(ee,2)>0);
    n2 = sum(sum(ee,1)>0);
    n_nodes = n1+n2;
    nedge = sum(sum(ee));
else
    n1 = sum(sum(ee,2)>0);
    n_nodes = n1;
    nedge = sum(sum(ee));
end

disp(['# of nodes: ' num2str(n_nodes)]);
disp(['# of edges: ' num2str(nedge)]);
end