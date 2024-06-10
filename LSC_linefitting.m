function [ label_Result1] = LSC_linefitting(data,numberOfModel,inlierScale,model_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[fitfn,resfn,degenfn,psize,numpar] = getModelPara(model_type);

kdtreeX = vl_kdtreebuild(data);
[neighborX, ~] = vl_kdtreequery(kdtreeX, data, data, 'NumNeighbors', psize) ;
N=size(data,2);
m=1;
res=[];inx=[];
for i=1:N
    pinx=neighborX(:,i);
    psub = data(:,pinx);
    % Check for degeneracy.
    isdegen = feval(degenfn,psub);
    if ~isdegen
        st = feval(fitfn,psub);
        % Compute residuals.
        ds = feval(resfn,st,data);
        inx(:,m)=pinx;
        res(:,m) = ds;
        m=m+1;
    end
end
R=res;

[N,M]=size(R);

inlier_index=double(R<inlierScale);
%---- step 1. compute preference matrix.
P = exp(-R/(inlierScale));

P=P.*inlier_index;
tic;
[U,S,V]=svds(P,numberOfModel);
toc
point=(U*S)'; % point index in the low dimension
point2=(V*S)';


weight=arrayfun(@(x)norm(point(:,x)),1:size(point,2));  % the norm of coordinate of each data point
weight2=arrayfun(@(x)norm(point2(:,x)),1:size(point2,2));
% figure;bar(weight);
% figure;hist(weight,50);
[II, EE]=entropy_thresholding(weight);
[II2, EE2]=entropy_thresholding(weight2);
keepInx_point_GMM=find(II>EE);
keepInx_point_GMM2=find(II2>EE2);
%[ keepInx_point_GMM, ~, ~, ~ ] = GMMremove(weight,GMMThreshold);
inlier=point(:,keepInx_point_GMM);
N2=size(inlier,2);
inlier2=point2(:,keepInx_point_GMM2);
R2=R(keepInx_point_GMM,:);
R3=R2(:,keepInx_point_GMM2);

inlier_index_2=double(R3<inlierScale);
%---- step 1. compute preference matrix.
P_2 = exp(-R3/(inlierScale));

P_2=P_2.*inlier_index_2;

%---- step 2. singular value decomposition of the preference matrix.
[U_2,S_2,V_2]=svds(P_2,numberOfModel);
point_2=(U_2*S_2)'; % point index in the low dimension

kdtreeX = vl_kdtreebuild(point_2);
[neighborX, ~] = vl_kdtreequery(kdtreeX, point_2, point_2, 'NumNeighbors', psize+2) ;
AA='xx';
[score,res_3,index_inx,scale] = hypothesis_updata(data,keepInx_point_GMM,neighborX,AA,model_type);

R_3=res_3;
[N_3,M_3]=size(R_3);

%sum(label_hypo_3>0)/M_3
inlier_index_3=double(R_3<inlierScale);
%---- step 1. compute preference matrix.
P_3 = exp(-R_3/(inlierScale));

P_3=P_3.*inlier_index_3;

%---- step 2. singular value decomposition of the preference matrix.
[U_3,S_3,V_3]=svds(P_3,numberOfModel);
point_3=(U_3*S_3)'; % point index in the low dimension
point2_3=(V_3*S_3)';


weight_3=arrayfun(@(x)norm(point_3(:,x)),1:size(point_3,2));  % the norm of coordinate of each data point
% figure;bar(weight);
% figure;hist(weight,50);
[II_3, EE_3]=entropy_thresholding2(weight_3);
keepInx_point_GMM_3=find(II_3>EE_3);

weight_3_2=arrayfun(@(x)norm(point2_3(:,x)),1:size(point2_3,2));
[II_3_2, EE_3_2]=entropy_thresholding2(weight_3_2);
keepInx_point_GMM_3_2=find(II_3_2>EE_3_2);

%%%%remove the hytpotheses with outliers
outliers=ones(1,N);
outliers(keepInx_point_GMM_3)=0;
current_inlier=find(outliers);
aa=sum(index_inx(current_inlier,:),1);
remove_label=find(aa>0);
[remove_hyp,ia,ib]=intersect(keepInx_point_GMM_3_2,remove_label);
keepInx_point_GMM_3_2(ia)=[];
%  [remove_hyp,ia,ib]=intersect(keepInx_point_GMM_3_2,remove_label);
%remove_hyp=find(keepInx_point_GMM_3_2==remove_label);
hypo_3_2=point2_3(:,keepInx_point_GMM_3_2);


ratioinliers=1;
% nbpoints=sum(keepInx_point_GMM); %11 13 1.5 9 8 3 5   9 13
[totm,totd]= ransac_sampling(hypo_3_2);
P = zeros(size(totd));
P(totd<ratioinliers) = 1;
[Q, ~] = refit_consensus_for_clsa(hypo_3_2, P, totm, ratioinliers);
F = pruning_consensus_for_clsa(Q);
%numberofmodel=1;
%select the k structures covering more points
W = ransacov_for_clsa(F,numberOfModel);
elabel=zeros(M_3,1);
for i = 1:size(W,2)
    elabel((W(:,i)>0))=i;
end
weight2=score;totdd=100*ones(N,10);outliers2=ones(1,N);
for i=1:numberOfModel
    remain_label=find(weight2>0);
    current_index=keepInx_point_GMM_3_2(find(elabel==i));
    [~,best_index]=max(weight2(current_index));
    best_index=current_index(best_index);
    rold=res_3(:,best_index);
    % Scales=newScales(best_index);
    newScales=scale(best_index);
    idinliers2=rold<newScales*2.5;
    idinliers=rold<newScales;
    outliers2(idinliers2)=0;
    aa=sum(index_inx(idinliers,:),1);

    totdd(:,i)=rold;

    % aa=sum(index_inx(current_inlier,:),1);
    remove_label=find(aa>0);
    remove_label=[remove_label current_index];
    weight2(remove_label)=-1;

end

[~,elabel]=min(totdd,[],2);
elabel(outliers2>0)=0;
elabel(outliers>0)=0;
label_Result1=elabel;
end

