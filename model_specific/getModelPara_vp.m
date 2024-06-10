function [ fitfn resfn degenfn psize numpar ] = getModelPara(model_type)

%---------------------------
% Model specific parameters.
%---------------------------

switch model_type
    
    case 'line'
        fitfn = @line_fit;
        resfn = @line_res;
        degenfn = @line_degen;
        psize = 2;
        numpar = 3;
     case '3Dline'
        fitfn = @line3D_fit;
        resfn = @line3D_res;
        degenfn = @line_degen;
        psize = 2;
        numpar = 6;
    case 'circle'
        fitfn = @fittingfn_circle2d;
        resfn = @distfn_circle2d;
        degenfn = @degenfn_circle2d;
        psize = 3;
        numpar = 3;
    case 'homography'
        fitfn = @homography_fit;
        resfn = @homography_res;
        degenfn = @homography_degen;
        psize = 4;
        numpar = 9;
    case 'fundamental'
        fitfn = @fundamental_fit;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 7;
        numpar = 9;
    case 'fundamental8'
        fitfn = @fundamental_fit8;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 8;
        numpar = 9;
    case 'motion'
        fitfn = @motion_fit;
        resfn = @motion_res;
        degenfn = @motion_degen;
        psize = 4;
        numpar = 25;
    case 'plane'
        fitfn=@plane_fit;
        resfn=@plane_res;
        degenfn=@plane_degen;
        psize=3;
        numpar=4;
    case 'vp'        
        fitfn = @fit_vp;
        resfn = @distSegVp;
        degenfn = @dummydegenerate;
        psize = 2;
        numpar=3;
    otherwise
        error('unknown model type!');
end

end