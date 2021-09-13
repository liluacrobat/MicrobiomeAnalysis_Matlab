function OTU_tbl = loadOTUtbl(tbl_name,tax_flag)
% read OTU table and extract sample id, taxonomy, and counts
if nargin<2
    tax_flag = 0;
else
    tax_flag = 1;
end
tbl = readtable(tbl_name,'delimiter','\t');
OTU_tbl.tax = table2array(tbl(:,1));
OTU_tbl.sample_id = tbl.Properties.VariableNames(2:end-tax_flag);
OTU_tbl.counts = table2array(tbl(:,2:end-tax_flag));
if tax_flag==1
    OTU_tbl.otu = OTU_tbl.tax;
    OTU_tbl.tax = table2array(tbl(:,end));
end
end