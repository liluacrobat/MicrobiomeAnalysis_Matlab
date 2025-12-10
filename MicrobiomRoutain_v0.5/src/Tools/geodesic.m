function D=geodesic(Data,K,distance,verbose,Para)
N = size(Data,1); 
if nargin<2
    K=1;
end
if nargin<3
    distance='cityblock';
end
if nargin<4
    verbose=1;
end
D=pdist2(Data,Data,distance);
%%%%% Step 1: Construct neighborhood graph %%%%%
disp('Constructing neighborhood graph...'); 
[~, ind] = sort(D); 
for i=1:N
  D(i,ind((2+K):end,i)) = inf; 
end
D = min(D,D');    %% Make sure distance matrix is symmetric
% Finite entries in D now correspond to distances between neighboring points. 
% Infinite entries (really, equal to INF) in D now correspond to 
%   non-neighoring points. 
%%%%% Step 2: Compute shortest paths %%%%%
disp('Computing shortest paths...'); 
% We use Floyd's algorithm, which produces the best performance in Matlab. 
% Dijkstra's algorithm is significantly more efficient for sparse graphs, 
% but requires for-loops that are very slow to run in Matlab.  A significantly 
% faster implementation of Isomap that calls a MEX file for Dijkstra's 
% algorithm can be found in isomap2.m (and the accompanying files
% dijkstra.c and dijkstra.dll). 
tic; 
for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
     if ((verbose == 1) && (rem(k,20) == 0)) 
          disp([' Iteration: ' num2str(k) '     Estimated time to completion: ' num2str((N-k)*toc/k/60) ' minutes']); 
     end
end
D(isinf(D))=0;
end