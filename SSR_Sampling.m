function [ label_Result1] = SSR_Sampling(data,numberOfModel,inlierScale,label,model_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figuredebug=0;

[fitfn,resfn,degenfn,psize,numpar] = getModelPara(model_type);

ikk=500;LSK=10;kernel_flag=1;
[~,R,inx] =random_sampling2(data,fitfn,resfn,degenfn,psize,numpar,ikk);
[N,M]=size(R);
label_hypo=zeros(1,ikk);
for i=1:ikk
    temp=label(inx(:,i));
    if(max(temp)==min(temp))
        label_hypo(i)= max(temp);
    end
end

inlier_index=double(R<inlierScale);
%---- step 1. compute preference matrix.
P = exp(-R/(inlierScale));

P=P.*inlier_index;

%---- step 2. singular value decomposition of the preference matrix.
[U,S,V]=svds(P,numberOfModel);

% [est_Phi, est_Theta, logPw_z]=learn_GibbsLDA(double(P), numberOfModel, inlierScale, inlierScale, 100, 10, 10, 1, [4 4]);
% Y = tsne(P,'Algorithm','exact','Standardize',true,'NumDimensions',numberOfModel,'Perplexity',20);
%rank(P)
%---- step 3. remove outliers.
% % % Sdiag= diag(S);
% % % Sdiag=1./Sdiag;
% % % S(logical(eye(numberOfModel)))=Sdiag;
point=(U*S)'; % point index in the low dimension
point2=(V*S)';
%point=Y';
%point=(P*V)';
if figuredebug
    colo=[0 0 1; 0 1 0; 1 0 1; 0 1 1; 1 0 0; 1 1 0; 0 0 0;  1 0.5 0];
    Youcolo=zeros(N,3);
    for inx_co=1:N
        Youcolo(inx_co,:)=colo(label(inx_co)+1,:);
    end
    
    if numberOfModel<3
        figure;gscatter(point(1,:),point(2,:),label);
    else
        figure;scatter3(point(1,:),point(2,:),point(3,:),35,Youcolo,'filled');%'b.'
    end
    for inx_co=1:M
        Youcolo2(inx_co,:)=colo(label_hypo(inx_co)+1,:);
    end
    figure;scatter3(point2(1,:),point2(2,:),point2(3,:),35,Youcolo2,'filled');%'b.'
end

% figure;scatter3(V(1,:),V(2,:),V(3,:),10,Youcolo2);
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
label2=label(keepInx_point_GMM);
N2=size(inlier,2);
if figuredebug
    Youcolo=zeros(N2,3);
    
    for inx_co=1:N2
        Youcolo(inx_co,:)=colo(label2(inx_co)+1,:);
    end
    
    if numberOfModel<3
        figure;gscatter(inlier(1,:),inlier(2,:),label2);
    else
        figure;scatter3(inlier(1,:),inlier(2,:),inlier(3,:),35,Youcolo,'filled');
    end
end
inlier2=point2(:,keepInx_point_GMM2);
label_hypo2=label_hypo(keepInx_point_GMM2);
M2=size(inlier2,2);
if figuredebug
    Youcolo=zeros(M2,3);
    
    for inx_co=1:M2
        Youcolo(inx_co,:)=colo(label_hypo2(inx_co)+1,:);
    end
    
    if numberOfModel<3
        figure;gscatter(inlier2(1,:),inlier2(2,:),label_hypo2);
    else
        figure;scatter3(inlier2(1,:),inlier2(2,:),inlier2(3,:),35,Youcolo,'filled');
    end
end
R2=R(keepInx_point_GMM,:);
R3=R2(:,keepInx_point_GMM2);

inlier_index_2=double(R3<inlierScale);
%---- step 1. compute preference matrix.
P_2 = exp(-R3/(inlierScale));

P_2=P_2.*inlier_index_2;

%---- step 2. singular value decomposition of the preference matrix.
[U_2,S_2,V_2]=svds(P_2,numberOfModel);
point_2=(U_2*S_2)'; % point index in the low dimension
point2_2=(V_2*S_2)';
if figuredebug
    colo=[0 0 1; 0 1 0; 1 0 1; 0 1 1; 1 0 0; 1 1 0; 0 0 0;  1 0.5 0];
    Youcolo=zeros(N2,3);
    for inx_co=1:N2
        Youcolo(inx_co,:)=colo(label2(inx_co)+1,:);
    end
    
    if numberOfModel<3
        figure;gscatter(point(1,:),point(2,:),label2);
    else
        figure;scatter3(point_2(1,:),point_2(2,:),point_2(3,:),35,Youcolo,'filled');%'b.'
    end
    Youcolo2=zeros(M2,3);
    for inx_co=1:M2
        Youcolo2(inx_co,:)=colo(label_hypo2(inx_co)+1,:);
    end
    figure;scatter3(point2_2(1,:),point2_2(2,:),point2_2(3,:),35,Youcolo2,'filled');%'b.'
end

% figure;
% plot(data(1,keepInx_point_GMM),data(2,keepInx_point_GMM),'b+');
%  axis ([0 10 0 10]);

%%%guided sampling

kdtreeX = vl_kdtreebuild(point_2);
[neighborX, ~] = vl_kdtreequery(kdtreeX, point_2, point_2, 'NumNeighbors', psize) ;
m=1;
param=[];res_3=[];inx_3=[];score=[];scale=[];index_inx=[];inlier_index=[];
h0=(243*0.6/(35*0.2*0.2*N2))^0.2;
for i=1:N2
    
    pinx=keepInx_point_GMM(neighborX(:,i));
    psub = data(:,pinx);
    % Check for degeneracy.
    isdegen = feval(degenfn,psub);
    if ~isdegen
        st = feval(fitfn,psub);
        % Compute residuals.
        ds = feval(resfn,st,data);
        st = reshape(st, numpar, 1);
        res_3(:,m) = ds;
        inx_3(:,m) = pinx;
        sub_index2=zeros(N,1);
        sub_index2(pinx)=1;
        index_inx(:,m)=sub_index2;
        %sub_index3=zeros(N,1);
        inlier_index(:,m)=ds<0.1;
        param(:,m)=st;
        [SRes,I] = sort(ds);
        scales_js=Hz_ILKOSE_NDF(SRes, LSK);
        delta=scales_js(end);
        hh=h0*delta;
        scale(m)=delta;
        score(m)=Hz_density_at_points_mkernel(SRes, 0, hh, kernel_flag);
        m=m+1;
    end
end
R_3=res_3;
[N_3,M_3]=size(R_3);

for ci=1:M_3
    residuals_Ci=res_3(:,ci);
    C=(scale(ci)).^2;
    Data_Prob_Peak_Ci=gaussian(residuals_Ci,0,C)'; %%%%The probability that data belonging to the ci_th sorted peaks.
    Data_Prob_Peaks(ci,:)=Data_Prob_Peak_Ci;
end


label_hypo_3=zeros(1,M_3);
for i=1:M_3
    temp=label(inx_3(:,i));
    if(max(temp)==min(temp))
        label_hypo_3(i)= max(temp);
    end
end
sum(label_hypo_3>0)/M_3
inlier_index_3=double(R_3<inlierScale);
%---- step 1. compute preference matrix.
P_3 = exp(-R_3/(inlierScale));

P_3=P_3.*inlier_index_3;

%---- step 2. singular value decomposition of the preference matrix.
[U_3,S_3,V_3]=svds(P_3,numberOfModel);
point_3=(U_3*S_3)'; % point index in the low dimension
point2_3=(V_3*S_3)';
if figuredebug
    colo=[0 0 1; 0 1 0; 1 0 1; 0 1 1; 1 0 0; 1 1 0; 0 0 0;  1 0.5 0];
    Youcolo_3=zeros(N,3);
    for inx_co=1:N
        Youcolo_3(inx_co,:)=colo(label(inx_co)+1,:);
    end
    
    if numberOfModel<3
        figure;gscatter(point_3(1,:),point_3(2,:),label);
    else
        figure;scatter3(point_3(1,:),point_3(2,:),point_3(3,:),35,Youcolo_3,'filled');%'b.'
    end
    Youcolo2_3=zeros(M_3,3);
    for inx_co=1:M_3
        Youcolo2_3(inx_co,:)=colo(label_hypo_3(inx_co)+1,:);
    end
    figure;scatter3(point2_3(1,:),point2_3(2,:),point2_3(3,:),35,Youcolo2_3,'filled');%'b.'
end

weight_3=arrayfun(@(x)norm(point_3(:,x)),1:size(point_3,2));  % the norm of coordinate of each data point
% figure;bar(weight);
% figure;hist(weight,50);
[II_3, EE_3]=entropy_thresholding2(weight_3);
keepInx_point_GMM_3=find(II_3>EE_3);
%[ keepInx_point_GMM, ~, ~, ~ ] = GMMremove(weight,GMMThreshold);
inlier_3=point(:,keepInx_point_GMM_3);
label2_3=label(keepInx_point_GMM_3);

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
label2_3_2=label_hypo_3(keepInx_point_GMM_3_2);

N2_3=size(inlier_3,2);
M_3_2=size(hypo_3_2,2);
if figuredebug
    Youcolo_3=zeros(N2_3,3);
    
    for inx_co=1:N2_3
        Youcolo_3(inx_co,:)=colo(label2_3(inx_co)+1,:);
    end
    
    if numberOfModel<3
        figure;gscatter(inlier_3(1,:),inlier_3(2,:),label2_3);
    else
        figure;scatter3(inlier_3(1,:),inlier_3(2,:),inlier_3(3,:),35,Youcolo_3,'filled');
    end
    Youcolo2_3=zeros(M_3_2,3);
    for inx_co=1:M_3_2
        Youcolo2_3(inx_co,:)=colo(label2_3_2(inx_co)+1,:);
    end
    figure;scatter3(hypo_3_2(1,:),hypo_3_2(2,:),hypo_3_2(3,:),35,Youcolo2_3,'filled');%'b.'
end



ratioinliers=1.5;
% nbpoints=sum(keepInx_point_GMM);
numberofSample=200;
tic;
[totm,totd]= ransac_sampling2(hypo_3_2,numberofSample);
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
if figuredebug
    Youcolo2_3=zeros(M_3_2,3);
    for inx_co=1:M_3_2
        Youcolo2_3(inx_co,:)=colo(elabel(inx_co)+1,:);
    end
    figure;scatter3(hypo_3_2(1,:),hypo_3_2(2,:),hypo_3_2(3,:),35,Youcolo2_3,'filled');%'b.'
end

% figure;
% plot(data(1,keepInx_point_GMM_3),data(2,keepInx_point_GMM_3),'b+');
%  axis ([0 10 0 10]);
%
%      figure;
%     for ij=1:M_3
%
%     line_plot(param(:,(ij)));hold on
%     axis ([0 10 0 10]);
%     end
%
%       figure;
%     for ij=1:numel(keepInx_point_GMM_3_2)
%
%     line_plot(param(:,keepInx_point_GMM_3_2(ij)));hold on
%     axis ([0 10 0 10]);
%     end

unique_elabel=unique(elabel);
estimated_hypo=[];elabel2=zeros(size(data,2),1);
weight2=score;totdd=100*ones(N,10);outliers=ones(1,N);
for i=1:numberOfModel
    remain_label=find(weight2>0);
    current_index=keepInx_point_GMM_3_2(find(elabel==i));
    %       figure;
    %     for ij=1:numel(current_index)
    %
    %     line_plot(param(:,current_index(ij)));hold on
    %     axis ([0 10 0 10]);
    %     end
    
    %current_index=remain_label;
    [value_,best_index]=max(weight2(current_index));
    best_index=current_index(best_index);
    %    estimated_hypo(:,i)=param(:,best_index);
    %     figure;
    %     line_plot(estimated_hypo(:,i));
    %     axis ([0 10 0 10]);
    rold=res_3(:,best_index);
    % Scales=newScales(best_index);
    current_inlier=rold<0.15;
    outliers(current_inlier)=0;
    totdd(:,i)=rold;
    
    %     n1=numel(remain_label);
    %         mu_matrix=ones(n1,1);
    %     Prob_data_Peak1=Data_Prob_Peaks(best_index,:);
    %     for ci=1:n1
    %         ci_label=remain_label(ci);
    %         Prob_data_Peak2=Data_Prob_Peaks(ci_label,:);
    %         mu_matrix(ci)=log(N*sum(Prob_data_Peak1.*Prob_data_Peak2)/(sum(Prob_data_Peak1)*sum(Prob_data_Peak2)));
    %     end
    %     aa=mu_matrix>0;
    
    
    
    % aa=sum(inlier_index(current_inlier,:),1);
    aa=sum(index_inx(current_inlier,:),1);
    remove_label=find(aa>0);
    %aa=mu_matrix(b,:)>0;
    %   figure;
    %     for ij=1:numel(remove_label)
    %
    %     line_plot(param(:,remove_label(ij)));hold on
    %     axis ([0 10 0 10]);
    %     end
    %cu_label=remain_label(aa>0);
    weight2(remove_label)=-1;
    
end

[~,elabel]=min(totdd,[],2);
elabel(outliers>0)=0;

% figure
% set(gca,'Position',[.02 .02 .95 .95]);
% gscatter(data(1,:),data(2,:),elabel,'bgmcryk','..........',20);
% axis ([0 10 0 10]);
% legend('off');
% hold on
% box off
% rr=axis;plot(rr(1:2),[rr(4), rr(4)],'k', [rr(2), rr(2)], rr(3:4),'k');
% set(gca,'XTick',[],'YTick',[]);
% box off
% ax2 = axes('Position',get(gca,'Position'),...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none',...
%     'XColor','k','YColor','k');
% set(ax2,'YTick', []);
% set(ax2,'XTick', []);
% box on
% hold off
%error=segmentationError(label,double(elabel))*100


% ratioinliers=1.5;
% % nbpoints=sum(keepInx_point_GMM);
% numberofSample=200;
% tic;
% [totm,totd]= ransac_sampling2(inlier_3,numberofSample);
% P = zeros(size(totd));
% P(totd<ratioinliers) = 1;
% [Q, ~] = refit_consensus_for_clsa(inlier_3, P, totm, ratioinliers);
% F = pruning_consensus_for_clsa(Q);
% %numberofmodel=1;
% %select the k structures covering more points
% W = ransacov_for_clsa(F,numberOfModel);
% elabel=zeros(N2_3,1);
% for i = 1:size(W,2)
%     elabel((W(:,i)>0))=i;
% end
% label_Result1=zeros(1,N);
% label_Result1(keepInx_point_GMM_3)=elabel;%label_inlier;
%
% error=segmentationError(label,double(label_Result1))*100
%
%


label_Result1=elabel;
% label_Result2=zeros(1,N);
% label_Result2(keepInx_point_GMM)=label_inlier;
end

