function [H,rms] = vgg_H_from_x_nonlin(H_initial,p1,p2)

% [H,rms] = vgg_H_from_x_nonlin(H_initial,xs1,xs2)
%
% Compute H using non-linear method which minimizes Sampson's approx to
% geometric reprojection error (see Hartley & Zisserman Alg 3.3 page 98 in
%                              1st edition, Alg 4.3 page 114 in 2nd edition). 

% An initial estimate of
% H is required, which would usually be obtained using
% vgg_H_from_x_linear. It is not necessary to precondition the
% supplied points.
%
% The format of the xs is
% [x1 x2 x3 ... xn ; 
%  y1 y2 y3 ... yn ;
%  w1 w2 w3 ... wn]

[r,c] = size(p1);

if (size(p1) ~= size(p2))
 error ('Input point sets are different sizes!')
end

global gp1 gp2 C1 C2;

gp1 = p1;
gp2 = p2;

% Make conditioners
C1 = vgg_conditioner_from_pts(p1);
C2 = vgg_conditioner_from_pts(p2);

H_cond = C2 * H_initial * inv(C1);

% opt = optimset( optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt','LargeScale','off', 'Diagnostics','off', 'Display','off','TolFun','1e-3');
opt = optimset('Diagnostics','off','Display','off','TolFun',1e-3);


% H_cond = lsqnonlin(@lsq_func,H_cond,[],[],opt);
H_cond = lsqnonlin(@lsq_func,H_cond,[],[],opt);

H = inv(C2) * H_cond * C1;
rms = sqrt(vgg_H_sampson_distance_sqr(H,gp1,gp2));

clear gp1 gp2 C1 C2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = lsq_func(H)
  global gp1 gp2 C1 C2;
  
  H_decond = inv(C2) * H * C1;
  
  d = vgg_H_sampson_distance_sqr(H_decond,gp1,gp2);
  
  r = sqrt(d);
  
