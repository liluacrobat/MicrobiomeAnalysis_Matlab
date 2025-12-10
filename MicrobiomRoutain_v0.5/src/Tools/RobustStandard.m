function [F,S]=RobustStandard(D,p,q)
if nargin<2
    p=0.1;
    q=1;
end
[d,N]=size(D);
M=mean(D,2);
D=D-repmat(M,1,N);
S=std(D,[],2);
SP=quantile(S(S~=0),p,1);
AS=S+SP*q;
F=D./repmat(AS,1,N);
F(S==0,:)=0;
end