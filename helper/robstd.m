function  s = robstd(x, type)
%Peter J. Rousseeuw & Christophe Croux scale estimators 

n = length(x);

A = abs(bsxfun(@minus,x,x'));

switch type
    case 's' %Sn
        s = c(n) * 1.1926*lhmedian(lhmedian(A,'hi'),'lo');
        %s = 1.1926*median(median(A))
        
    case 'q' %Qn
        y = A(triu(true(size(A)),1));
        
        s = d(n) * 2.2219 * quantile(y,.25);
        
    case 'm' %MAD
        s = 1.4826*mad(x,1);
        
    otherwise
        error('unrecognized option');
        
end


end


%----- auxiliary function

function  c = c(n)
% finite sample correction  for Sn

table = [1 0.743 1.851  0.954  1.351  0.993 1.198 1.005 1.131];

if logical(rem(n,2))
    
    c = n/(n-0.9); %dispari
else
    c = 1; %pari
end

if n < 9
    c = table(n);
end
end


function  d = d(n)
% finite sample correction  for Qn

table = [1 0.399 0.994  0.512  0.844  0.611 0.857 0.669 0.872];

if logical(rem(n,2))
    
    d= n/(n+1.4); %dispsri
else
    d = n/(n+3.8); %pari
end

if n < 9
    d = table(n);
end
end



function  m = lhmedian(A,type)
% low and high median used in Sn


n = length(A);

B = sort(A);

switch type
    case 'lo'
        m = B(floor((n+1)/2));
    case 'hi'
        m = B(floor(n/2)+1,:);
    otherwise
        error('unrecognized option');
end
end

