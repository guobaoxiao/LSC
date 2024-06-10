function [totd,totm, mss, tim ] = mssWeighted_for_all( X, m, blksz, distance, model, win, epsi )
%MSSWEIGHTED Summary of this function goes here
%   Detailed explanation goes here

if (~exist('win','var') || isempty(win))
    win = 0.1;
end

if(~exist('epsi','var') || isempty(epsi))
    tau=1;
end

[distFun, hpFun, fit_model, cardmss, isdegen, d] = set_model( model );

if (strcmp(model,'subspace4'))
    d = size(X,1)*4;
elseif (strcmp(model,'subspace9'))
    d = size(X,1)*9;
elseif (strcmp(model,'affspace3'))
    d=size(X,1)*4;
end


n = size(X,2);
mss = nan(m,cardmss);
H_cap = nan(d,m);
R_cap = nan(n,m);

tic

for i = 1:m
    if(i <= blksz)
        w = ones(1,n); % initialize weight
        seedinx = mod(i,n)+1;   % seed
        mss(i,1) = seedinx;     % store
        w(seedinx) = 0;         % no replacement
        % populate i-th mss
        j = 2;
        while(j<=cardmss)
            othinx = randsample(n,1,true,w);     % extract a point
            mss(i,j) = othinx;                   % store
            w(othinx) = 0;                       % no replacement
            %------ ceck for degeneracy ------%           
            if(~isdegen( X(:,mss(i,1:j)) ) )
                j = j+1;
            elseif(all(w==0))
                error('I can not find a non degenerate mss');    
            end
            %---------------------------------%      
        end  
    else
        if (mod(i-1,blksz)==0)
            %use the previous block in order to update sampling wheights
            k = (i-1)/blksz;
            new_blk = (k-1)*blksz+1:i-1;
            H_cap(:,new_blk) = hpFun( X, mss(new_blk,:) );
            R_cap(:,new_blk) = res( X, H_cap(:,new_blk), distFun );
            switch distance
                
                case 'uni'
                    W =  ones(n,n);
                    
                case 'cos'
                    D = squareform(pdist(R_cap(:,1:i-1),'cosine'));
                    sigma = quantile(D(:),win)/3;
                    W =  exp(-(D).^2./sigma^2);
                    
                case 'exp'
                    D = squareform(pdist(exp(-R_cap(:,1:i-1)),@tanimoto));
                    sigma = quantile(D(:),win)/3;
                    W =  exp(-(D).^2./sigma^2);
                    
                case 'loc'
                    D = squareform(pdist(X'));
                    sigma = quantile(D(:),win)/3;
                    W =  exp(-(D).^2./sigma^2);
                case 'mgs'
                    [~, resinx]  = sort(R_cap(:,1:i-1),2);
                    W = computeIntersection(resinx',resinx',win*blksz);     
                case 'tani'
                    tau=epsi/5;
                    D = squareform(pdist(exp(-R_cap(:,1:i-1)./tau),@tanimoto));
                    sigma = quantile(D(:),win)/3;
                    W =  exp(-(D).^2./sigma^2);
                    
                case 'cauchy'
                    
                    D = squareform(pdist(prefMat(R_cap(:,1:i-1),epsi, 6),@tanimoto));
                    sigma = quantile(D(:),win)/3;
                    W =  exp(-(D).^2./sigma^2);
                    
            end
            
        end
        
        
        seedinx = mod(i,n)+1;   % seed
        w = W(seedinx,:);       % initialize weight
        mss(i,1) = seedinx;     % store
        w(seedinx) = 0;         % no replacement
        
        
        % populate i-th mss
        j = 2;
        while(j<=cardmss)
            
            if (all(w==0))
                othinx = randsample(n,1);         % extract a point
            else
                othinx = randsample(n,1,true,w);  % extract a point
            end
            
            mss(i,j) = othinx;                   % store
            w(othinx) = 0;                       % no replacement
            w = w.*w;                            % update weight  
            %--------- ceck for degeneracy ---------%   
            if( ~isdegen( X(:,mss(i,1:j)) ) )
                j = j+1;
            elseif(all(w==0))
                disp('I can not find a non degenerate mss');
            end
            %---------------------------------------%           
        end      
    end  
end
mss=mss(blksz+1:end,:);
totm = hpFun(X,mss); %hypotheses
totd = Dres( X, totm, distFun ); disp('Residuals computed');
tim = toc;
%fprintf(strcat(distance,': elapsed time for sampling %i\n'),tim);









