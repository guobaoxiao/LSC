function  C = seg_from_U( U,X, model, out)
%SEG_FROM_H Summary of this function goes here
%   Detailed explanation goes here

if(nargin<4)
    out=1; %oulier are present in the data
end

[ distFun, hpFun, fit_model, cardmss, idegen, d] = set_model( model );
N = size(U,1); % number of points
k = size(U,2); % number of models


[mss, C] = mss_from_U( U , cardmss);
M = hpFun(X,mss); % models
R = res( X, M, distFun );

if (out)
    % compute the inlier threshold
    epsilon = nan(k,1);
    theta = 3.5;
    
    for i = 1: k
        r= R(C==i,i);
        epsilon(i) =  x84thresh( r, theta );
        
        current_outliers = C==i & R(:,i)>epsilon(i);
        C(current_outliers) = 0 ;
        
    end
    
    
MEend


end




