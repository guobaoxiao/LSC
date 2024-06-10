function [inx,res]= ransac_sampling(data)
% Prepare storage.
%-----------------
[dis,npts] = size(data);
% par = zeros(numpar,M);
 res = zeros(npts,npts);
 inx = zeros(dis,npts);

    
for  m=1:npts

        pinx = m;%randsample(npts,1);
        psub = data(:,pinx);

    ds = feval(@higdis_line_res,psub,data);

%     par(:,m) = st;
    res(:,m) = ds;
    inx(:,m) = psub;
end






