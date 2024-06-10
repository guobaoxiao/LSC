function [II, EE]=entropy_thresholding(Ws)
%%%%%From "A Density-Based Data Reduction Algorithm for Robust Estimators"
%%%%%Lecture Notes In Computer Science; Vol. 4478, 2007

%%%%%Input a set of weights of particles; 
%%%%%Output the threshold EE and the values used for the threshold II>EE.
%qq0=Ws.^2; 
qq0=Ws; 
%qq=max(qq0)-qq0; 
%qq=median(qq0)-qq0; 
qq=mean(qq0)-qq0; 
qq(qq<0)=0;
%qq=exp(-qq0);
%  qq=sqrt(max(qq0)-qq0); 
  %qq=(max(qq0)-qq0).^2; 
 %qq=qq0;
pp=qq/sum(qq); 
II=-log(pp+eps); 
EE=sum(pp.*II); 
% big_weight_Index=find(II>EE)

% ABs_BigWeights3=[As(big_weight_Index)' Bs(big_weight_Index)' Ws(big_weight_Index)' Scales(big_weight_Index)']';
% figure(200);clf; plot3(ABs_BigWeights3(1,:), ABs_BigWeights3(2,:), ABs_BigWeights3(3,:)./ABs_BigWeights3(4,:), 'r*');  axis(xy_range); view(93,8); xlabel('A', 'FontWeight', 'bold', 'FontSize', 16); ylabel('B', 'FontWeight', 'bold', 'FontSize', 16); zlabel('Weight', 'FontWeight', 'bold', 'FontSize', 16); 
% figure(201); clf; plot(ABs_BigWeights3(1,:), ABs_BigWeights3(2,:),'r*'); xlabel('A', 'FontWeight', 'bold', 'FontSize', 16); ylabel('B', 'FontWeight', 'bold', 'FontSize', 16); %axis(xy_range)
% 
% exp(-II)