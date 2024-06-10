function C = segmentation( mss, X, model, U , sigma, cost ,str)
%SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here




[ distFun, hpFun, fit_model, cardmss, idegen, d] = set_model( model );


n = size(X,2); % number of points
k = size(mss,1); % number of models



% initialization
% (1) fit models on estimated mss

M = hpFun(X, mss); % models
R = Dres( X, M, distFun );


%  (2) initialize segmentation

if (strcmp(str, 'nearest'))
    C = nan(n,1);
    for i = 1 : n
        [~, C(i)]=min( R(i,:) );
    end
    
elseif(strcmp(str, 'fromu'))
    C = seg_from_binaryU(U);
end

%  make sure at least cardmss inliers are assigned
for i=1:k
    C(mss(i,:))=i;
end


% (3) inlier threshold & rejects outliers


tune = 5;



epsilon = nan(k,1);
for i = 1: k
    
    r= R(C==i,i);
    r=r(r<5*sigma);
    epsilon(i) = cost * tune *  robstd(r, 's');
    C((R(:,i) > epsilon(i))& C==i)=0;
end
%figure; imagesc(C); pause;

% (4) refit models

Mls = fit_model(X(:,C>0), C(C>0));
Rls = Dres( X, Mls, distFun );

% (5) inlier threshold and outlier rejection
epsilon = nan(k,1);
for i = 1: k
    r= Rls(C==i,i);
    r=r(r<5*sigma);
    epsilon(i) = cost * tune *  robstd(r, 's');
    C((R(:,i) > epsilon(i))& C==i)=0;
end


% (6) deal with intersection
for i = 1 : n
    [t, u ] = min( Rls(i,:) );
    if(t <= epsilon(u))
        C(i)=u;
    else
        C(i)=0;
    end
end


%figure; imagesc(C)
end







