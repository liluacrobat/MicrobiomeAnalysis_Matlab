function dis=BioDis(X,Y,method)
switch method
    case 'BC'
        dis=BC_dis(X,Y);
    case 'KL'
        dist=KLDiv(X,Y);
    case 'Hellinger'
        dis=Hellinger_dis(X,Y);
    case 'BRAYCD'
        dis=braycd(X,Y);
end
end
function dis=BC_dis(X,Y)
X=X/sum(X);
Y=Y/sum(Y);
BC=sum(sqrt(X.*Y));
dis=-log(BC);
end
function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end
% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
end
end
function dis=Hellinger_dis(X,Y)
X=X/sum(X);
Y=Y/sum(Y);
BC=sum(sqrt(X.*Y));
dis=sqrt(1-BC);
end
function D = braycd(X1,X2)
% created by Olivier Zemb
% input : X1 and X2 as two row vectors containing the observations of the
% site 1 and site 2
% output : the bray curtis distance
% Note that this function can be used in combination with pdist.
% braycd= (Si+Sj -2 Cij)/(Si+Sj)
X1=X1/sum(X1);
X2=X2/sum(X2);
X=[X1; X2];
Cij=sum(min(X));
Si_plus_Sj=sum(X1)+sum(X2);
D=(Si_plus_Sj-2*Cij)/Si_plus_Sj;
end


