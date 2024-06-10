close all;
clear;
clc;

model_type = 'circle';
[fitfn,resfn,degenfn,psize,numpar] = getModelPara(model_type);


%A='Pami4_threeCircle';
%A='Pami4_fourCircle';
%A='Pami4_fiveCircle';
A='Pami4_sixCircle';
load(['data\', A,'.mat']);
%save(['D:\data\', A,'.mat'],'data','label');

groundtrue=label;
sigma=1.5;
GMMthreshold=1;

t=0.1;


tic
X=data;
numberOfModel=numel(unique(groundtrue))-1;
label=groundtrue';
[ labelResult1] = LSC_linefitting(data,numberOfModel,sigma,model_type);
time=toc;
elabel=zeros(size(data,2),1);
%inlier_index=ones(size(data,2),1);

for i_model=1:numberOfModel
    current_index=find(labelResult1==i_model);
    %sample_index=unidrnd(numel(current_index),1,psize);
    param=feval(fitfn,data(:,current_index));
    rold= feval(resfn, param, data);
    %     sr=sort(abs(rold));
    %     scales_js=Hz_ILKOSE_NDF(sr, LS_K);
    %     delta=scales_js(end)
    elabel(rold<t)=i_model;
    %inlier_index(rold<1.0)=0;
    %final_res(:,i_model)=rold;
end


error=segmentationError(label',double(elabel'))*100;
fprintf('t=%.2f;sigma=%.2f:%s: mean std time %.3f %.3f %.3f \n',t,sigma,A,mean(error),std(error),mean(time));
elabel2=elabel;
for i=1:numberOfModel
    temp=find(elabel==i);
    temp2=tabulate(label(temp));
    [a,b]=max(temp2(:,2));
    elabel2(temp)=temp2(b,1);
end

colo=[0 0 1; 0 1 1;0 1 0; 1 0 1;  1 0 0; 1 1 0; 0 0 0;  1 0.5 0];
figure
set(gca,'Position',[.02 .02 .95 .95]);
gscatter(X(1,:),X(2,:),label,colo,'............',25);
%axis ([0 100 0 100]);
legend('off');
hold on
box off
rr=axis;plot(rr(1:2),[rr(4), rr(4)],'k', [rr(2), rr(2)], rr(3:4),'k');
set(gca,'XTick',[],'YTick',[]);
box off
ax2 = axes('Position',get(gca,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
hold off

figure
set(gca,'Position',[.02 .02 .95 .95]);
gscatter(X(1,:),X(2,:),elabel2,colo,'............',25);
%axis ([0 100 0 100]);
legend('off');
hold on
box off
rr=axis;plot(rr(1:2),[rr(4), rr(4)],'k', [rr(2), rr(2)], rr(3:4),'k');
set(gca,'XTick',[],'YTick',[]);
box off
ax2 = axes('Position',get(gca,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
hold off
% figure;
% bar(weight);
%
% figure;
% gscatter(density,delta,label(labelInlier),colo,'x.....',10)








