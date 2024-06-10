function [ X0, me ] = guess_init( X, k ,g )
%GUESS_INIT guess initialization for nnmf


MAXiter = 2000; % Maximum number of iterations for KMeans 
REPlic = 600; % Number of replications for KMeans
c  = kmeans(X,k,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
X0  = clust_indicator( c );

if nargin==3
    me = Misclassification(c(g>0),g(g>0));
end


end

