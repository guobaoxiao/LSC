%------------------------------------------------------------------------
% Function to calculate distances between a 2d line and an array of points.
function [dist, P] = higdis_line_res(P, X)  
    X=X';
    P=P';
    

% mindis = bsxfun(@minus, X, P);
%     dist=sqrt(sum(mindis.^2,2));
%     
%     
    cos=X*P'./(sqrt(sum(X.^2,2)).*sqrt(sum(P.^2,2)));
    dist=sqrt(sum(X.^2,2)).*(1.-cos.^2);
    dist(dist<0)=0;
end