%======================================================================

clear all;
close all;

data_path = 'E:\Hgeneratedata';%eneratedata
data_files = dir(data_path);
data_files(1:2) = [];
num=length(data_files);
%test_label=find(Final_H_results_SE>10);

for ii = 1:num%numel(test_label)
    %     data = char(dataList(i));
    %ii=test_label(iij);
    %fprintf('--------------%d-----------\n',ii);
    %model_type = 'fundamental8';
    model_type = 'homography';
    debug=0;
    Dataname=data_files(ii).name;
    file_name = strcat(data_path, '/', Dataname);
    AA=Dataname(1:end-4);
    pack=load(file_name);
    X = pack.data;
    groundtrue=pack.label;
    scores=pack.score;
    I1 = pack.img1;
    I2 = pack.img2;

    [x1 T1] = normalise2dpts(X(1:3,:));
    [x2 T2] = normalise2dpts(X(4:6,:));
    data = [ x1 ; x2 ];

    ratioinliers=0.1;

    sigma=0.05;
    tic;
    %elabel = SSR_single_homogr(data,sigma,LSK,k,model_type);
    elabel = LSC_SH(data,sigma,ratioinliers,model_type);
    ttime=toc;

    error1=segmentationError(groundtrue,(elabel'))*100;
    RFM_SCAN_results(ii,:)=[error1,ttime];
    fprintf('SSR: %s error:%f;  Time: %f\n', AA,error1,ttime);

end
fprintf('LSC: error:%f; Time: %f\n',mean(RFM_SCAN_results(:,1)),mean(RFM_SCAN_results(:,2)));


