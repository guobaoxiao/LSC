function P = line_fit(X)
N = size(X, 2);
A = [X; ones(1, N)];
[XYn, T] = normalise2dpts(A);


[u d v] = svd(XYn',0);   % Singular value decomposition.
P = v(:,3);              % Solution is last column of v.

% Denormalise the solution
P = T'*P;


%[U,S,V] = svd(A);
%P = V(:, 3);
end