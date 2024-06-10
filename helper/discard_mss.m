function [ impure ] = discard_mss( S, F )
%DISCARD_MSS Summary discard possibly impure mss according to segmentation F:
% impure is a boolean vector: true if the corresponding mss is impure
%  A mss is considered pure if:
%         - all its points have the same F-labels~=0;
%         -  in case of more labels there is a predominant one

[m, cardmss] = size(S);
quorum = ceil(cardmss/2)+1;
impure = false(m,1);




for i= 1:m
    labels = unique(F(S(i,:)));
    num_labels = numel(labels);
    
    if (num_labels>1)
        
        % counting the number of points with labels(i)
        frq = nan(1, num_labels); % frequency
        for j = 1:num_labels
            frq(j)= sum(F(S(i,:))==labels(j));
        end
        
        if(max(frq)<quorum)
            impure(i) = true;
        end
        
    end
    
end

% figure;
% subplot(1,4,1); imagesc(F(S)); title('mss')
% subplot(1,4,2); imagesc(impure); title('impure')
% subplot(1,4,3); imagesc(F(S(impure,:))); title('mss pure')
% subplot(1,4,4); imagesc(F(S(impure,:))); title('mss impure')

end

