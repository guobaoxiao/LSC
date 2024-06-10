function [ Q ,S ] = rinforzino(X, S, P, F, softIndU , model, sigma, niter, G)
%RINFORZINO exploit clustering information contained in F: replace impure samples with pure
%ones



[ distFun, hpFun, fit_model, cardmss, idegen, d] = set_model( model );


% disacard impute mss


% [~, F_hp] = countMss( S, F );
% impure = F_hp==0;
impure  = discard_mss( S, F );
Q = P(:, ~impure);
S = S(~impure,:);
%Q = P;




Z = mss_from_softU( softIndU, cardmss, niter );
S = [S;Z];
H = hpFun(X, Z); % models
R = Dres( X, H, distFun );
Q  = [Q, prefMat(R, sigma, 6)];





if(nargin>8)
    
    [f1, ~] = countMss( S, G );
    [f2 , ~] = countMss( [S;Z], G );
    figure;
    subplot(1,2,1); bar(f1); title('#mss per models'); subplot(1,2,2); bar(f2); title('#mss per models')
end

