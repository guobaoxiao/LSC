function  I  = clust_indicator( C )
%CLUST_INDICATOR nxk matrix n number of points k number of cluster;
% (i,j)=1 if point i belongs to cluster j
n = numel(C);
k = max(C);
I = zeros(n,k);
for i = 1:k
    I(C == i,i) = 1;
    
end
end

