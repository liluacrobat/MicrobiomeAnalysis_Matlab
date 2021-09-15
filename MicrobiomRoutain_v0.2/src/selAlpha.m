function y = selAlpha(x,sel)

y = x;
if isempty(x.ob_otu)
    y.ob_otu = [];
else
    y.ob_otu = x.ob_otu(sel);
end

if isempty(x.shannon)
    y.shannon = [];
else
    y.shannon = x.shannon(sel);
end

if isempty(x.chao1)
    y.chao1 = [];
else
    y.chao1 = x.chao1(sel);
end
end
