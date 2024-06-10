function P = flat_fit(A)

% Fit 2D affine subspace.
m = size(A,1);
n = size(A,2);
mu = sum(A,2)./n;
Ab = A - repmat(mu,1,n);
[ U ~, ~, ] = svd(Ab);
B = U(:,1:2);
P = reshape([mu B],m*3,1);
    
end

