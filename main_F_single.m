%======================================================================
clear all;
close all;
data_path = 'E:\Fgeneratedata';
data_files = dir(data_path);
data_files(1:2) = [];

n=numel(data_files);


RFM_SCAN_results=[];
for ii = 1:n
    %     data = char(dataList(i));
    %close all;
    % fprintf('--------------%d-----------\n',ii);
    model_type = 'fundamental8';
    %model_type = 'homography';
    debug=0;
    Dataname=data_files(ii).name;
    file_name = strcat(data_path, '/', Dataname);
    AA=Dataname(1:end-4);
    pack=load(file_name);
    X = pack.data;
    I1 = pack.img1;
    I2 = pack.img2;
    xy=X;
    groundtrue=pack.label;
    scores=pack.score;
    I1 = pack.img1;
    I2 = pack.img2;

    [x1 T1] = normalise2dpts(X(1:3,:));
    [x2 T2] = normalise2dpts(X(4:6,:));
    data = [ x1 ; x2 ];

    nbpoints = size(X(1:3,:), 2);
    ratioinliers=0.01;
    [fitfn,resfn,degenfn,psize,numpar] = getModelPara(model_type);
    tic;
    elabel = LSC_SH(data,ratioinliers,ratioinliers,model_type);

    ttime=toc;
    error1=segmentationError(groundtrue,(elabel'))*100;
    fprintf('LSC: ii=%d  error:%f;  Time: %f\n',ii,error1,ttime);

end



