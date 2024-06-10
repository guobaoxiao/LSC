function H = hpVp( X, S )
%HPVP Summary of this function goes here
%   Detailed explanation goes here
m = size(S,1);
H = zeros(3,m);
for i = 1 : m
    
    
    % data to fit
    x1 = X(:,S(i,1));
    x2 = X(:,S(i,2));
    l1 =  cross([x1([1,3]);1], [x1([2,4]);1]);
    l2 = cross([x2([1,3]);1], [x2([2,4]);1]);
    a =  cross(l1,l2);
    if(all(a==0))
        a=rand(3,1);
    end
    H(:,i)=    a;%./a(3);
    
    
    
    %     figure
    %     line(x1([1,2]),x1([3,4]),'linewidth',5);
    %     hold on
    %     line(x2([1,2]),x2([3,4]),'linewidth',5);
    %     pause
    %     hold on;
    %     plot(H(1,i),H(2,i), 'r*','MarkerSize',10)
    
    
    
end

end

