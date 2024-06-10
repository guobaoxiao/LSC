function [ Q, K ] = refit_consensus_for_clsa( X, P, H, epsilon )

%REFIT_PREF refit preferences, accept a refit only if the CS is increased
%   at the moment work only for hard preferences
%   INPUT:
%          X            data points
%          epsilon      inlier threshold
%          P            preference matrix
%          H            sampled hypotheses
%   OUTPUT:
%          Q            new preference matrix
%          K            new sampled hypotheses

%[fit_model,distFun,degenfn,cardmss,numpar] = getModelPara(model_type);
m   = size(P,2); % number of hypotheses

count = 0;
if (numel(epsilon)==1)
    epsilon=epsilon*ones(1,m);
end
Q = P;
K = H;
for j = 1:m
    % compute consensus set
    inlier  = P(:,j)>0;
    card_cs = sum(inlier); % cardinality of the consensus set
    % refit model
    if(card_cs>1)
        h_new  = mean(X(:, inlier),2); % new hypothesis
        % update if necessary
         r=feval(@higdis_line_res, h_new, X);
        %r = res(X, h_new, distFun);
        if(sum(r<epsilon(j)) >= card_cs)
            K(:,j) = h_new;
            Q(:,j) = prefMat(r, epsilon(j), 0 );
            count  = count+1;
        end
        
    end
end

%fprintf('Number of updated hypotheses %i\n', count);
end

