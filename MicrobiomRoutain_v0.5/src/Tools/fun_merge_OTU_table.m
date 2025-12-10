function R = fun_merge_OTU_table(x)
tax = [];
otu = [];
for i=1:length(x)
    otu=[otu;x{i}.otu];
    tax = [tax;x{i}.tax];
    [otu,ia,~] = unique(otu);
    tax = tax(ia);
end
[u_otu,ia,~] = unique(otu);
u_tax = tax(ia);
mtx = zeros(length(u_otu),length(x));
for i=1:length(x)
    ia = AlignMem(u_otu,x{i}.otu);
    mtx(ia~=0,i) = x{i}.mtx(ia(ia~=0),1);
end
R.mtx = mtx;
R.otu = u_otu;
R.tax = u_tax;
end