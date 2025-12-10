function adjp = fun_bonferroni(p)
%% adjust for multiple testing correction uisng Bonferroni method
p_vertor = p(:);
adjp = p_vertor*length(p_vertor);
adjp(adjp>1) = 1; 
adjp = reshape(adjp, size(p));
end