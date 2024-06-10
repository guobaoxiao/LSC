function  P  = prefMat( R, epsilon, r, tau )
%
% prefMat Estimate the preference matrix frome the residual matrix using
% a voting procedure specified by r.
% 
% INPUT:
%       R: residual matrix: the entry (i,j) is the residual between the i-th
%       point and the j-th model
%       epsilon: threshold value
%       r: specifies the voting scheme:
%           r=1 exponential voting
%           r=2 gaussian voting
%           r=0 hard voting
%           r=3 tukey
%           r=4 nocutoff
%
%OUTPUT:
%       P: preference matrix, each row depicts points preferences with
%       repsect to the putative model hypotheses
%
% Author: Luca Magri
% For any comment, question or suggestion about the code please contact
% luca (dot) magri (at) unimi (dot) it


[n,m]=size(R);
%initialize Preference Matrix
P = zeros(n,m);
%set cutoff value
% if ~exist('r','var') && ~exist('tau','var')
%     r=1;
%     tau=epsilon/5;
%     %disp('Exponential voting is used')
% end

switch r
    case 0
        %disp('Hard voting')
        I=R<epsilon;
        P(I) = 1;
    case 1
        disp('Exponential voting')
        %tau=epsilon/5;
        I=R<epsilon;
        P(I) = exp(-R(I)./tau);
    case 2
        disp('Gaussian voting')
        sigma=epsilon/4;
        I=R<epsilon;
        P(I) = exp(-(R(I).^2)./(sigma^2));
    case 3
        disp('Tukey voting')
        I=R<epsilon;
        P(I) = 1-tukey(R(I),epsilon);
    case 4
        disp('No cutoff');
        %tau= 2.5*mad(R(:),1)/5;
        tau= 3.5*std(R(:))/5;
        P = exp(-R./tau);
    case 5
        %disp('Time costant')
        tau = epsilon/5;
        P = exp(-R./tau);
    case 6
        c = 5* epsilon; % epsilon is the std
        P = 1./(1+(R./c).^2);
        
end




end

