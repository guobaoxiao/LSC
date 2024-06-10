function [ U ] = pruning_consensus( P,cardmss )
%pruning_hedge pruning hyperedge
%  extract from the preference matrix hyperdeges.
% P,is interpreted as the incidence matrix of an hypergraph
% i)   largest hyperdeges
% ii)  not contained in the union of the selected hyperedge (rationale: at least a new point is explained)
% iii) covering the nodes (possibly with redundancy)



if nargin==1
    cardmss=0;
end

% step 0:  initialization
[n,m]   = size(P);
row_flg = zeros(n,1);
col_flg = ones(1,m);
col_idx = 1:m;
U       = [];

% step 1: ordering columns
sumP = sum(P,1);
[~, order] = sort(sumP, 'descend');
P = P(:, order);


% step 2: main loop
union = zeros(n,1);

while(any(row_flg==0) && any(col_flg==1)) % until all the points are covered
    
    % keep the largest hyperedge among the available ones
    col_available = col_idx(col_flg==1);
    pick = col_available(1); %index of the first available column
    col_flg(pick) = 0;
    U = [U, P(:, pick)];
    k = size(U,2);
    row_flg(U(:,k)==1)=1; % deflagga i punti
    % delete intersection
    j = k;
    union = sum(U(:,1:k),2);
    union(union>0)=1;
    while(j<m)
        j = j+1;
        if(col_flg(j))
            % se è un sottoinsieme di U(:,k) deflagga
            % se è contenuto nell'unione
            if(all (P(:,j).* union== P(:,j)))
                col_flg(j) = 0;
            end
        end
    end
    
    
end

U=U(:,sum(U,1)>cardmss);

r=(sum(U,2)==0);

for i=1:numel(r)
    if(r(i)==1)
        t=zeros(n,1);
        t(i)=1;
        U=[U,t];
    end
end


end

