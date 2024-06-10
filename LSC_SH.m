function [ label_Result1] = LSC_SH(data,inlierScale,ratioinliers,model_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[fitfn,resfn,degenfn,psize,numpar] = getModelPara(model_type);
%LSK=10;
%%%%%%%%%%%ininitsampling
kdtreeX = vl_kdtreebuild(data);
[neighborX, ~] = vl_kdtreequery(kdtreeX, data, data, 'NumNeighbors', psize) ;
N=size(data,2);
m=1;
res=[];
parfor i=1:N
    pinx=neighborX(:,i);
    psub = data(:,pinx);
    % Check for degeneracy.
    isdegen = feval(degenfn,psub);
    if ~isdegen
        st = feval(fitfn,psub);
        % Compute residuals.
        ds = feval(resfn,st,data);
        res(:,i) = ds;
        %m=m+1;
    end
end
res(:,sum(res,2)==0)=[];
R=res;

% if m==1
%     label_Result1=zeros(N,1);
%     return;
% end

inlier_index=double(R<inlierScale);
%---- step 1. compute preference matrix.
P = exp(-R/(inlierScale));

P=P.*inlier_index;

%---- step 2. singular value decomposition of the preference matrix.
[U,S,V]=svds(P,1);
point=(U*S)'; % point index in the low dimension

% figure;scatter3(V(1,:),V(2,:),V(3,:),10,Youcolo2);
weight=arrayfun(@(x)norm(point(:,x)),1:size(point,2));  % the norm of coordinate of each data point

[II, EE]=entropy_thresholding(weight);

keepInx_point_GMM=find(II>EE);

%[ keepInx_point_GMM, ~, ~, ~ ] = GMMremove(weight,GMMThreshold);
inlier=point(:,keepInx_point_GMM);

point_2=inlier;
kdtreeX = vl_kdtreebuild(point_2);
[neighborX, ~] = vl_kdtreequery(kdtreeX, point_2, point_2, 'NumNeighbors', psize+2) ;
AA='hs';
[score,res_3,index_inx,scale] = hypothesis_updata(data,keepInx_point_GMM,neighborX,AA,model_type);

[value_,best_index]=max(score);
R_3=res_3(:,best_index);
elabel=R_3<ratioinliers;


label_Result1=elabel;% label_Result2=zeros(1,N);

% label_Result2(keepInx_point_GMM)=label_inlier;
end

