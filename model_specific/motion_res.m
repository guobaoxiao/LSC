%%%motion segmentation -abao-2014-11-17

function dist = motion_res(D,indx)

%dist=sum((P*x).^2);
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
    U = V(:, 1:c); 
 D = bsxfun(@minus, D, mean(D_sub)); 
newX = (eye(F) - U*U')*D'; 
dist = sum(newX.^2, 1); 
% U=U(:,1:c);
%         P=squeeze(U);
%         K=size(P,1);
%         P=eye(K)-P*P';
%         dist=sum((P*D').^2);
% f = size(X,1);
% U = reshape(U,f,4);
% % proietto  X su U ottenendo Y
% k=X'*U; %coefficienti di Y in U
% Y = k(1).*U(:,1)+k(2).*U(:,2)+k(3).*U(:,3)+k(4).*U(:,4);
% dist= norm(X-Y);
end


    