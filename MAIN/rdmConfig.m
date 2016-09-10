function [config] = rdmConfig(target, template, Vtarget,  Vtemplate, ... 
    motion, normspec, anisotropy, scale, lambda, ctrlptflag)
%RDML2 will match template onto target through a transformation of type
%'motion'

if(~exist('ctrlptflag', 'var'))
    ctrlptflag = 1;
end

config.template = template;
%handle empty case
config.Vtemplate = Vtemplate;

config.target = target;
%handle empty case
config.Vtarget = Vtarget;

%handle empty case
config.motion = motion;
config.normspec = normspec;
config.specstring = strcat(motion, '-', anisotropy, '-', normspec);

[n1,d] = size(template);
[n2,d] = size(target);

config.nTemp = n1;
config.nTarg = n2;
config.scale = scale;
config.lambda = lambda;

config.display = 0;
config.init_param = [ ];
config.max_iter = 500;
config.normalize = 0;

if ctrlptflag == 1
    if(d==3)
        config.ctrl_pts = createCtrlPts(template, .5, .5, .5, 6);
    else
        config.ctrl_pts = createCtrlPts(template, .15, .15, 5); 
    end
else
    config.ctrl_pts = createCtrlPts(template,1);
end

config.alpha = 1;

config.beta = 1/n1^2;0;
config.opt_affine = 1;

[n,d] = size(config.ctrl_pts); 

init_affine = repmat([zeros(1,d) 1],1,d);

config.init_affine = init_affine;

config.tempWts = ones(1, n1)/n1;
config.targWts = ones(1, n2)/n2;


config.gradflag = 1;

config.Lb_affine = -inf*repmat([ones(1,d) 1],1,d);
config.Ub_affine = inf*repmat([ones(1,d) 1],1,d);
config.Lb_tps = -inf*ones(n-d-1,d);
config.Ub_tps = inf*ones(n-d-1,d);

config.Lb_flips = -ones(1,n2);
config.Ub_flips = ones(1,n2);


config.Ub_pa = ones(1,n2);
config.Lb_pa = zeros(1,n2);

switch lower(config.motion)
    case 'tps'

    config.init_nonrigid = zeros(n-d-1,d);
    config.oldP =  zeros(d, (n-d-1))';

        switch lower(config.normspec)
            case 'restricted'
                config.init_flips = .1*ones(1,n2);
                config.init_param = [init_affine reshape(config.init_nonrigid, 1, d*n-d*(d+1)) config.init_flips];

            case 'fixed'
                config.init_param = [init_affine reshape(config.init_nonrigid, 1, d*n-d*(d+1))];
                
            case 'free'
                config.init_param = [init_affine reshape(config.init_nonrigid, 1, d*n-d*(d+1)) 0*config.Vtarget'];
                %do it this way for the normals: forget about the 
                %parameters and just add the result to Vtarget.
        end
        
    case 'grbf'
        config.init_nonrigid = zeros(n,d);
        config.oldP =  zeros(d, n)';
        
        switch lower(config.normspec)
            case 'restricted'
                config.init_flips = .1*ones(1,n2);
                config.init_param = [reshape(config.init_nonrigid, 1, d*n) config.init_flips];
            case 'fixed'
                config.init_param = [reshape(config.init_nonrigid, 1, d*n)];
            case 'free'
                config.init_param = [reshape(config.init_nonrigid, 1, d*n) 0*config.Vtarget'];
        
            
        end
    otherwise
        [x0,Lb,Ub] = set_bounds(motion, n1);
        config.init_param = x0;
        config.Lb = Lb;
        config.Ub = Ub;
        config.gradflag = 1;
end

