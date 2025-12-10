function SampleSize_Estimation
% Hajian-Tilaki, K., 2014. Sample size estimation in diagnostic test studies 
% of biomedical informatics. Journal of biomedical informatics, 48, pp.193-204.
i=0;
for auc = 0.75:0.001:0.95
    i=i+1;
    Prev = 0.5;
    d=0.075;
    a = norminv(auc)*1.414;
    % n2 = n*Prev;
    % n1 = n*(1-Prev);
    nV = (0.0099*exp(-a^2/2))*((5*a^2+8)/Prev+(a^2+8)/(1-Prev));
    n(i) = 1.96^2*nV/d^2;
end
n=[];
i=0;
for d = 0.05:0.0005:0.2
    i=i+1;
    Prev = 0.5;
    auc = 0.75;
    a = norminv(auc)*1.414;
    % n2 = n*Prev;
    % n1 = n*(1-Prev);
    nV = (0.0099*exp(-a^2/2))*((5*a^2+8)/Prev+(a^2+8)/(1-Prev));
    n(i) = 1.96^2*nV/d^2;
end
a=(0.05:0.0005:0.2)*2;
end