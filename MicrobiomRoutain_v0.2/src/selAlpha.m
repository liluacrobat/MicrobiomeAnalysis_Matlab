function y = selAlpha(x,sel)
y = x;
y.ob_otu = x.ob_otu(sel);
y.shannon = x.shannon(sel);
y.simposon = x.simposon(sel);
y.chao1 = x.chao1(sel);
end