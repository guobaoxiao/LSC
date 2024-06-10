function display_model_vp(model_type,file,X,labels,parameters,groundtruth,img1,img2)
%f
%
X = X';
[ fitfn,resfn,degenfn,psize,numpar] = getModelPara_vp(model_type);
switch model_type
    case 'line'
        figure(3);
        subplot(1,2,1);
        plot(X(1,:),X(2,:),'k.','MarkerSize',10);
        gscatter(X(1,:), X(2,:), labels', 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        title('The first label in red is outlier label');
        
        subplot(1,2,2);
        gscatter(X(1,:), X(2,:), labels', 'rgbcmykgbcmk', '+osdv*<>phosd', 10);
        title('The first label in red is outlier label');
        
        figure(4);
        gscatter(X(1,:), X(2,:), labels', 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);hold on;
        box on;
        
        for i = 1:size(parameters,1)
            line_plot(parameters(i,:),i);hold on;
            %             axis([-1,1,-1,1]);
        end
        
    case 'circle'
        figure(3);
        subplot(1,2,1);
        % imshow(img1);hold on;
        plot(X(1,:),X(2,:),'k.','MarkerSize',10);
        % gscatter(X(1,:), X(2,:), AKSWP_T2, 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        
        subplot(1,2,2);
        % imshow(img2);hold on;
        gscatter(X(1,:), X(2,:), labels, 'rgbcmykgbcmk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        
        figure(4);
        % gscatter(X(1,:), X(2,:), AKSWP_T2, 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);hold on;
        box on;
        
        for i = 1:size(parameters,1)
            visualfn_circle2d(parameters(i,1:numpar));hold on;
            %     axis([-1,1,-1,1]);
        end
        
    case 'homography'
        figure(3);
        subplot(1,2,1);
        imshow(img1);hold on;
        gscatter(X(1,:), X(2,:), labels, 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        
        subplot(1,2,2);
        imshow(img2);hold on;
        gscatter(X(4,:), X(5,:), labels, 'rgbcmykgbcmk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        if nargin == 7
            figure(4);
            subplot(1,2,1);
            title('Groundtruth img1');
            imshow(img1);hold on;
            gscatter(X(1,:),X(2,:),groundtruth,'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
            
            subplot(1,2,2);
            title('Groundtruth img2');
            imshow(img2);hold on;
            gscatter(X(4,:),X(5,:),groundtruth,'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        end
    case 'fundamental'
        figure(3);
        subplot(1,2,1);
        imshow(img1);hold on;
        gscatter(X(1,:), X(2,:), labels, 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        
        subplot(1,2,2);
        imshow(img2);hold on;
        gscatter(X(4,:), X(5,:), labels, 'rgbcmykgbcmk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        if nargin == 7
            figure(4);
            subplot(1,2,1);
            title('Groundtruth img1');
            imshow(img1);hold on;
            gscatter(X(1,:),X(2,:),groundtruth,'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
            
            subplot(1,2,2);
            title('Groundtruth img2');
            imshow(img2);hold on;
            gscatter(X(4,:),X(5,:),groundtruth,'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        end
    case 'fundamental8'
        figure(3);
        subplot(1,2,1);
        imshow(img1);hold on;
        gscatter(X(1,:), X(2,:), labels, 'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        
        subplot(1,2,2);
        imshow(img2);hold on;
        gscatter(X(4,:), X(5,:), labels, 'rgbcmykgbcmk', '+osdv*<>phosd', 10);
        % title('The first label in red is outlier label');
        if nargin == 7
            figure(4);
            subplot(1,2,1);
            title('Groundtruth img1');
            imshow(img1);hold on;
            gscatter(X(1,:),X(2,:),groundtruth,'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
            
            subplot(1,2,2);
            title('Groundtruth img2');
            imshow(img2);hold on;
            gscatter(X(4,:),X(5,:),groundtruth,'rgbcmykgbcmyk', '+osdv*<>phosd', 10);
        end
    case 'vp'
        %%%%%%%%%%%%%%%%%%%%%%%%% display of vp %%%%%%%%%%%%%%%%%%%%%%%%%
        color = {'r','g','y','b','m','c'};
        data = X;
        figure('position',[300 200 845.5 600]);set(gca,'Position',[.056 .10 .90 .85]);
%         subplot(1,2,1);
        imshow(img1,'border','tight','initialmagnification','fit');hold on;
        for c = min(groundtruth):max(groundtruth)
            inx = find(groundtruth==c);
            %     c_lines = data(:,inx);
            for ic = 1:length(inx)
                plot(data(1:2,inx(ic)),data(3:4,inx(ic)),char(color(c)),'LineWidth',2);
            end
        end
%        saveas(gca,strcat('E:\code\T-linkage\Tlinkage_beta\YorkUrbanDB\Figure\',file,'-1.eps'),'epsc2');
 %       saveas(gca,strcat('E:\code\T-linkage\Tlinkage_beta\YorkUrbanDB\Figure\',file,'-1.png'),'png');
        figure('position',[300 200 845.5 600]);set(gca,'Position',[.056 .10 .90 .85]);
        %         subplot(1,2,2);
        imshow(img1,'border','tight','initialmagnification','fit');hold on;
        cls = unique(labels);
        for c = 1:length(cls)
            inx = find(labels==cls(c));
            %     c_lines = data(:,inx);
            for ic = 1:length(inx)
                plot(data(1:2,inx(ic)),data(3:4,inx(ic)),char(color(c)),'LineWidth',2);
            end
        end
%        saveas(gca,strcat('E:\code\T-linkage\Tlinkage_beta\YorkUrbanDB\Figure\',file,'-2.eps'),'epsc2');
%        saveas(gca,strcat('E:\code\T-linkage\Tlinkage_beta\YorkUrbanDB\Figure\',file,'-2.png'),'png');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end