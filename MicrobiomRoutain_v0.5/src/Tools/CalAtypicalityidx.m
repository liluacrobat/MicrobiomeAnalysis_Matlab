function J = CalAtypicalityidx(Xref,Xtar)
X = [Xref,Xtar];
Rel_ref = CalRel(Xref);
[~,n] = size(X);
s = sum(Xref>0,2);
idx = s==size(Xref,2) & mean(Rel_ref,2)>0.0002;

X = X(idx,:);

Ref = CLR(Xref(idx,:));
Target = CLR(Xtar(idx,:));
[m,N] = size(Ref);
n = N-1;
rng default
D = 5;
d = D-1;
J = zeros(1000,size(Target,2));
for k = 1:1000
    r = randperm(m,5);
    Ref_i = Ref(r,:);
    Target_i = Target(r,:);
    Theta = cov(Ref_i');
    Mu = mean(Ref_i,2);
    for i = 1:size(Target,2)
        X_i = Target_i(:,i);
        qy = N/(N+1)*(X_i-Mu)'/Theta*(X_i-Mu);
        J(k,i) = betainc(qy/(qy+n),d/2,(n-d+1)/2);
    end
end
end