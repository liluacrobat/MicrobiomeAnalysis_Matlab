function f_ls = addPrefix(id_ls,prefix,num_flag)
if nargin<3
    num_flag=0;
end

n = length(id_ls);
f_ls = cell(n,1);
if num_flag==1
    
    for i=1:n
        f_ls{i} = strcat(prefix,num2str(id_ls(i)));
    end
else
    for i=1:n
        f_ls{i} = strcat(prefix,id_ls(i));
    end
    
end
end