function [ S ] = mss_from_softU( U, cardmss, numhp )
%MSS_FROM_SOFTU sampling mss exploiting columns of U as wehight
%   numhp is the number of hypotheses per segment


[n,k] = size(U);
S = nan(numhp*k, cardmss);


count = 1;

for i = 1 : k
    for j = 1 : numhp
    w = U(:,i);
    %seedinx = mod (count, n) + 1;
    
    
    for l = 1 : cardmss
        inx = randsample(n, 1, true, w);
        S(count, l)= inx;
        w(inx)=0;
        w=w.*w;
    end
    
    count = count + 1;
    end
    
end

