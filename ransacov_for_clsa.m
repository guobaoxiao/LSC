function V = ransacov_for_clsa( F, kappa)
% RANSACOV 
% Input: F a characteristic matrix, whose columns represent a family of
% set.
% If F is the only input, the minimum number of subsets from F that covers 
% the data is returned via the SET COVER formulation
% If the number kappa of desired set is also specified, the MAXIMUM
% COVERAGE problem is solved.
% Output: V contains as column the indicator functions of the attained sets.   


% When using the code in your research work, please cite the following paper:
% Luca Magri, Andrea Fusiello,  Multiple Models Fitting as a Set Coverage Problem, CVPR, 2016.
% For any comments, questions or suggestions about the code please contact
% magrilucal (at) gmail (dot) com


N = size(F,1);
r = size(F,2);

if(nargin==1)
    disp('set cover')
    
    A = -F;
    b = -ones(N,1);
    lb = zeros(r,1);
    f = ones(1,r);
    intcon = [1:r];
    
    
    x = intlinprog(f,intcon, A,b,[],[],lb);
    
    
    V=F(:,x>0);
else
    %%  maximum coverage
    
        %disp('maximum coverage')
        w=ones(1,N);

    
    f = [zeros(1,r),-w];
    
    Aeq=[ones(1,r), zeros(1,N)];
    beq = kappa;
    A =[-F, eye(N)];
    b = zeros(N,1);
    lb = zeros(r+N,1);
    ub = ones(r+N, 1);
    
    intcon = [1:N+r];
    x = intlinprog(f ,intcon, A, b, Aeq, beq, lb, ub);
    
    x = x(1:r);
    V=F(:,x>0);
end
end

