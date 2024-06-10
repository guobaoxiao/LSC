function [ label_Result1] = LSC_F(data,numberOfModel,inlierScale,AA,model_type)
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
for i=1:N
    pinx=neighborX(:,i);
    psub = data(:,pinx);
    % Check for degeneracy.
    isdegen = feval(degenfn,psub);
    if ~isdegen
        st = feval(fitfn,psub);
        % Compute residuals.
        ds = feval(resfn,st,data);
        res(:,m) = ds;
        m=m+1;
    end
end
R=res;
inlier_index=double(R<inlierScale);
%---- step 1. compute preference matrix.
P = exp(-R/(inlierScale));

P=P.*inlier_index;

%---- step 2. singular value decomposition of the preference matrix.
[U,S,V]=svds(P,numberOfModel);

point=(U*S)'; % point index in the low dimension

% figure;scatter3(V(1,:),V(2,:),V(3,:),10,Youcolo2);
weight=arrayfun(@(x)norm(point(:,x)),1:size(point,2));  % the norm of coordinate of each data point

[II, EE]=entropy_thresholding(weight);

keepInx_point_GMM=find(II>EE);

inlier=point(:,keepInx_point_GMM);

point_2=inlier; % point index in the low dimension

kdtreeX = vl_kdtreebuild(point_2);
[neighborX, ~] = vl_kdtreequery(kdtreeX, point_2, point_2, 'NumNeighbors', psize+2) ;
[score,res_3,index_inx,scale] = hypothesis_updata(data,keepInx_point_GMM,neighborX,AA,model_type);

%[score,res_3,index_inx,inx,scale] = hypothesis_updata(data,inx_3,LSK,k,n_iterations,model_type);
R_3=res_3;

inlier_index_3=double(R_3<inlierScale);
%---- step 1. compute preference matrix.
P_3 = exp(-R_3/(inlierScale));

P_3=P_3.*inlier_index_3;

%---- step 2. singular value decomposition of the preference matrix.
[U_3,S_3,V_3]=svds(P_3,numberOfModel);
point_3=(U_3*S_3)'; % point index in the low dimension
point2_3=(V_3*S_3)';

weight_3=arrayfun(@(x)norm(point_3(:,x)),1:size(point_3,2));  % the norm of coordinate of each data point

[II_3, EE_3]=entropy_thresholding2(weight_3);
keepInx_point_GMM_3=find(II_3>EE_3);

weight_3_2=arrayfun(@(x)norm(point2_3(:,x)),1:size(point2_3,2));
[II_3_2, EE_3_2]=entropy_thresholding(weight_3_2);
keepInx_point_GMM_3_2=find(II_3_2>EE_3_2);

% %%%%remove the hytpotheses with outliers
outlierss_label=ones(1,N);
outlierss_label(keepInx_point_GMM_3)=0;
hypo_3_2=point2_3(:,keepInx_point_GMM_3_2);
M_3_2=size(hypo_3_2,2);

ratioinliers=0.8;
[totm,totd]= ransac_sampling(hypo_3_2);
P = zeros(size(totd));
P(totd<ratioinliers) = 1;
[Q, ~] = refit_consensus_for_clsa(hypo_3_2, P, totm, ratioinliers);
F = pruning_consensus_for_clsa(Q);
%numberofmodel=1;
%select the k structures covering more points
W = ransacov_for_clsa(F,numberOfModel);
elabel=zeros(M_3_2,1);
for i = 1:size(W,2)
    elabel((W(:,i)>0))=i;
end

weight2=score;totdd=100*ones(N,10);outliers=ones(1,N);
for i=1:numberOfModel
    %remain_label=find(weight2>0);
    current_index=keepInx_point_GMM_3_2(find(elabel==i));
    %current_index=remain_label;
    [~,best_index]=max(weight2(current_index));
    best_index=current_index(best_index);
    rold=res_3(:,best_index);
    % Scales=newScales(best_index);
    newScales=scale(best_index);
    if newScales>0.1
        newScales=0.01;
    end
    idinliers2=rold<newScales*2.5*100;
    idinliers=rold<newScales*2.5*10;
    outliers(idinliers2)=0;
    aa=sum(index_inx(idinliers,:),1);
    remove_label=find(aa>2);
    totdd(:,i)=rold;
    weight2(remove_label)=-1;
end

[~,elabel]=min(totdd,[],2);
elabel(outliers>0)=0;
elabel(outlierss_label>0)=0;
label_Result1=elabel;
end

