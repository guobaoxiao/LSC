%%%motion segmentation -abao-2014-11-17

function P = motion_fit(D,indx)

% y=gramsmithorth(x);
% P=eye(size(x,1))-y*y';
%     [U,~,~]=svd(x);
%     O=U(:,1:4);
%     %O = orth(X(:,S(i,:)));
%    P=O(:); 
    c = 5; 
    F = size(D, 2); 
    c = min(c, numel(indx));
    D_sub = D(indx, :); 
    X = bsxfun(@minus, D_sub, mean(D_sub)); 
    [U, S, V] = svd(X, 'econ'); 
    
    S = diag(S); 
    id = find(cumsum(S)/sum(S) < 0.995); 
    if numel(id) 
        if (id(end) > c); 
            c = id(end);
        end
    else
        c = 1; 
    end
    P = V(:, 1:c); 
end


    