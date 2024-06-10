function [ M ] = fit_vp( X, C )
%FIT_VP expoiting the dual space


if nargin==1 || isempty(C)
    C=ones(size(X,2),1);
end

label=sort(unique(C));
kappa=length(label); %numero di cluster;
M=nan(3,kappa);



for i=1:kappa
    
    L  = label(i);
    segments = X(:,C==L);
   A=[];
   B=[];
    for j = 1:size(segments, 2)
         a = cross([segments([1,3],j);1],[segments([2,4],j);1]);% line = point to fit in dual space
         if(abs(a(3))>1e-3);
         a = a./sqrt(a(1)^2+a(2)^2);
         A = [A; a(1),a(2)];
         B=[B;-a(3)];
         end
    end

    U = A\B;
    M(:,i)=[U;1];
    
end


end

