function [ C ] = seg_from_binaryU( U )
%SEG_FROM_BINARYU Summary of this function goes here
%   Detailed explanation goes here


k = size(U,2);  % number of points
n = size(U,1);  % number of models

C=nan(n,1);

h = U;

indU =indMax(U);
for i = 1 : k
    C(indU(:,i)==1)=i;
end

end




