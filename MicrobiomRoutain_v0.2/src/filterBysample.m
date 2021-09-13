function f_tbl = filterBysample(OTU_tbl,idx)
if length(unique(idx))~=2
    idx = idx(idx~=0);
end
f_tbl.sample_id = OTU_tbl.sample_id(idx);
counts = OTU_tbl.counts(:,idx);
otu_sel = sum(counts,2)>0;
f_tbl.counts = counts(otu_sel,:);
f_tbl.tax = OTU_tbl.tax(otu_sel);
if isfield(OTU_tbl,'otu')
    f_tbl.otu = OTU_tbl.otu(otu_sel);
end
if isfield(OTU_tbl,'clr')
    f_tbl.clr = OTU_tbl.clr(otu_sel,idx);
end
if isfield(OTU_tbl,'rel')
    f_tbl.rel = OTU_tbl.rel(otu_sel,idx);
end
end