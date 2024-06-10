function res = flat_res(P,X)

% Compute residual to 2D affine subspace.
m = length(P)/3;
P = reshape(P,m,3);
mu = P(:,1);
B = P(:,2:3);
n = size(X,2);
res = sqrt(sum((X - ((B*B')*(X - repmat(mu,1,n)) + repmat(mu,1,n))).^2))';

end