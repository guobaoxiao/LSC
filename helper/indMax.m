function [ B ] = indMax( A )
%INDMAX max sulle righe
B = zeros(size(A));
[~,I] = max(A, [], 2); 
B(sub2ind(size(A), 1:length(I), I')) = 1;
B=logical(B);
end

