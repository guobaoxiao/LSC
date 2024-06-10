function [X,G,img1] = load_data_FM( sequencename, out )


if (nargin<2)
    % load inlier data
    out = 0;
end

load(sequencename)
img1=img1;

  X = data; G = label'; N=size(X,2);
 if(out==0)
     disp('only inlier are loaded')
     X = X(:,G~=0);
     G = G(G~=0);
 end
 

[X,flg] = remove_repeated_points(X);
G=G(flg);

% reorder data
[~, order]=sort(G);
G=G(order);
X=X(:,order);

end

