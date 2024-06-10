clear;
close all;
clc;

%-------------------------------------------------------------------------%

%--Load data--------------------------------------------------------------%
data_path = 'data\AdelaideRMF\H';
%data_path = 'data\AdelaideRMF\';
data_files = dir(data_path);
data_files(1:2) = [];
model_type='homography';
%model_type='fundamental8';

savename=[];
%    DD=[];
%for i=1:11
lsc_results=[];
[fitfn,resfn,degenfn,psize,numpar] = getModelPara(model_type);

j=1;
for kkk1 = 1:numel(data_files)

    Dataname=data_files(kkk1).name;
    file_name = strcat(data_path, '/', Dataname);
    pack=load(file_name);
    xy = pack.data;
    groundtrue=pack.label;
    [dat_img_1 T1] = normalise2dpts(xy(1:3,:));
    [dat_img_2 T2] = normalise2dpts(xy(4:6,:));
    data = [ dat_img_1 ; dat_img_2 ];
    X=data(1:2,:);
    Y=data(4:5,:);
    numberOfModel=numel(unique(groundtrue))-1;
    AA=Dataname(1:end-4);
    %-------------------------------------------------------------------------%

    %--- Prepare data --------------------------------------------------------%
    xy = pack.data;
    I1 = pack.img1;
    I2 = pack.img2;
    %scores=pack.score;
    groundtrue=pack.label;
    %---Normalize data--------------------------------------------------------%
    [dat_img_1 T1] = normalise2dpts(xy(1:3,:));
    [dat_img_2 T2] = normalise2dpts(xy(4:6,:));
    data = [ dat_img_1 ; dat_img_2 ];
    label=groundtrue';

    sigma=0.05; %0.04
    %k=ikk+6;
    if AA=="bonhall" || AA=="ladysymon"
        sigma=0.1;
    end
    if AA=="johnsonb"%9.21 19
        sigma=0.09;
    end
    tic

    elabel=LSC_H(data,numberOfModel,sigma,AA,model_type);
    time1=toc;
    error1=segmentationError(groundtrue,double(elabel'))*100;
    lsc_results(j,:)=[error1,time1];
    fprintf('Done %s iteration.....%.2f  %.2f\n', AA,error1,time1);
    j=j+1;
    [XX,Y,Z]=size(I1);
end
ER_LSC=lsc_results(:,1);
time_LSC=lsc_results(:,2);

RFM_MP=mean(lsc_results)
RFM_MP=median(lsc_results)






