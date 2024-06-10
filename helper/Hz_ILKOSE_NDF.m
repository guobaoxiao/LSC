function scales_js=Hz_ILKOSE_NDF(sr, LS_K)
%%%%%cumulative density function.%%%%%
 K=ceil(length(sr)*LS_K/100); k=K/length(sr);
 residual_K=sr(K); G_threhsold=2.5; 
 if k<1
 scale0=residual_K/norminv((1+k)/2);
 else
     scale0=sqrt(sum(sr.^2)/length(sr));
 end
 scalestmp=scale0; scales_js(1)=scale0; loop_flag=1;
while loop_flag
 sr_j=sr(sr<scalestmp*G_threhsold);
 if length(sr_j)<=0
    % disp('length(sr_j)=0'); 
     break;
 end
 k=K/length(sr_j);
  if k<1
    scale_j=residual_K/norminv((1+k)/2);
    if abs(scale_j-scalestmp)<=eps
        loop_flag=0; 
    else
        scalestmp=scale_j;
        scales_js(end+1)=scale_j; 
    end
  else
        break
        disp(['debug while. number_used_ratio_j=' num2str(k)]);
  end
end






% number_used=ceil(length(sr)*LS_K/100); number_used_ratio=number_used/length(sr);
%  scales=sr(number_used)/norminv((1+number_used_ratio)/2);
% %%%%try to use loop to refine the biased scale.
% scalestmp=scales; scales_js(1)=scales; loop_flag=1; k=0; 
% while loop_flag
%     sr_j=sr(sr<scalestmp*2.5); inlier_number_j=length(sr_j);
%     number_used_j=ceil(inlier_number_j*LS_K/100);
%     if number_used_j<1
%         disp('debug');
%     end
%     number_used_ratio_j=number_used_j/inlier_number_j;
%     scales_j=sr_j(number_used_j)/(norminv((1+number_used_ratio_j)/2)+eps);
%     if abs(scales_j-scalestmp)<0.01|scales_j>scalestmp|k>=100
%         loop_flag=0;% scales_js(end+1)=scales_j; 
%     else
%         scalestmp=scales_j; 
%         scales_js(end+1)=scales_j; k=k+1 ;
%     end
% end
% if k==100; 
%     k, scales_js(1:100)
% end
%     
% k
% scales_js

