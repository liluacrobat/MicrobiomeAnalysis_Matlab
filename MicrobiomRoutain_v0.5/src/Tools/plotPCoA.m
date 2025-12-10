function [scores,dis] = plotPCoA(rel,Label,d,dis)
% dis    = symmetric dissimilarity matrix
%
% ----- Methods: -----
%  'ave': Average distance                                             (=  'd2')
%   'bc': Bray-Curtis dissimilarity
%  'bin': Binomial deviance dissimilarity (scale-invariant)
% 'binU': Binomial deviance dissimilarity (unscaled)
%  'can': Canberra metric                                              (= 'd10')
% 'can2': Canberra dissimilarity
%  'chi': Chi-square metric                                            (= 'd15')
% 'chi2': Chi-square distance                                          (= 'd16')
%  'chJ': Chao's abundance-based Jaccard dissimilarity
% 'chJB': Bias-corrected version of 'chJ'
%  'chS': Chao's abundance-based Sorensen dissimilarity
%  'cho': Orloci's Chord distance                                      (=  'd3')
%  'cod': Coefficient of Divergence + exclude double-zeros             (=  d11')
%  'cor': Pearson Correlation dissimilarity
% 'cor2': Spearman Rank Correlation dissimilarity
%  'cy':  Chao's dissimilarity for count data
% 'czek': Czekanowski distance + exclude double-zeros                  (=  'd8')
%  'euc': Euclidean distance (default)                                 (=  'd1')
% 'euc2': Euclidean distance + exclude double-zeros
%  'esq': Euclidean distance squared
% 'eucs': Euclidean similarity
%  'geo': Geodesic metric                                              (=  'd4')
%  'gow': Gower's metric                                               (= 's15')
% 'gow2': Gower's metric + exclude double-zero                         (= 's19')
% 'gow3': Gower's metric + exclude double-zeros (log2 transform)
% 'gow4': Gower's metric + exclude double-zeros (log10 transform)
%  'hel': Hellinger dissimilarity                                      (= 'd17')
%  'ioa': Whittaker's Index of Association dissimilarity               (=  'd9')
%  'jac': Jaccard dissimilarity                                        (=  's7')
%  'kul': Kulczynski quantitative dissimilarity                        (= 's18')
%  'man': Manhattan (City-block) metric                                (=  'd7')
% 'man2': Manhattan (City-block) metric + exclude double-zeros         (=  'd8')
% 'mink': Minkowski's metric                                           (=  'd6')
%  'mor': Morisita's Index of Overlap for count data
% 'morH': Morisita-Horn dissimilarity
%  'och': Ochiai quantitative dissimilarity
%  'sor': Sorensen's dissimilarity                                     (=  's8')
%  'wat': Watson's nonmetric coefficient dissimilarity                 (= 'd13')
%   's1': Simple Matching dissimilarity
%   's2': Rogers & Tanimoto dissimilarity
%   's3': Sokal & Sneath dissimilarity
%   's4': Sokal & Sneath dissimilarity
%   's5': Sokal & Sneath dissimilarity
%   's6': Sokal & Sneath dissimilarity
%   's9': Jaccard dissimilarity variant
%  's10': Sokal & Sneath dissimilarity
%  's11': Russel & Rao dissimilarity
%  's13': Binary Kulczynski dissimilarity
%  's14': Binary Ochiai dissimilarity
%  's26': Faith dissimilarity
%  'spe': Species Profiles distance
rel = rel(sum(rel,2)~=0,:);
if nargin<4
    dis='euc';
end
if strcmpi(dis,'yc')==1
    dis = f_dis_yc(rel);
else
    dis = f_dis(rel',dis);
end
pcoa = f_pcoa(dis,1);
scores = pcoa.scores(:,1:d)';
power = pcoa.expl(:,1);
plotFIGURE(pcoa.scores(:,1:d)',Label,defaultColor);
set(gca,'FontSize',14);
xlabel(['PCo1 (' num2str(round(power(1)*1000)/10) '%)']);
ylabel(['PCo2 (' num2str(round(power(2)*1000)/10) '%)']);

pbaspect([1 1 1]);
end