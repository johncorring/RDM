function [optStruct] = createOptStruct(config)

optStruct = struct();

template = config.template;

[n,d] = size(config.ctrl_pts);
[m,d] = size(config.template);
[m1,d] = size(config.target);


%get input OPSs
optStruct.template = config.template;
optStruct.target = config.target;

optStruct.Vtemplate = config.Vtemplate;
optStruct.Vtarget = config.Vtarget;

optStruct.tempWts = config.tempWts;
optStruct.targWts = config.targWts;

optStruct.ctrlPts = config.ctrl_pts;

%get params
optStruct.scale = config.scale;
optStruct.beta = config.beta;
optStruct.lambda = config.lambda;
optStruct.d = d;
optStruct.n = n;

optStruct.oldP = config.oldP;
optStruct.motion = config.motion;




if(strcmp(config.motion, 'grbf'))
    
    [K,U,dUx,dUy,dUz] = compute_kernel2(config.ctrl_pts, template, 0, config.motion);
    x0 = [reshape(config.init_nonrigid', 1, numel(config.init_nonrigid))];
    optStruct.kernel = K;
    optStruct.K = K;
    optStruct.U = U;
    optStruct.basis = U;
    optStruct.Tx = dUx;
    optStruct.Ty = dUy;
    optStruct.Tz = dUz;
elseif(strcmp(config.motion, 'tps'))
    [K,U,dUx,dUy,dUz] = compute_kernel2(config.ctrl_pts, template, 0, config.motion);
    
    Pm = [ones(m,1) template];
    Pn = [ones(n,1) config.ctrl_pts];
    PP = null(Pn');  % or use qr(Pn)
    basis = [Pm U*PP];
    
    if(d == 2)
        Tx = [zeros(m,1) ones(m,1) zeros(m,1) dUx*PP];
        Ty = [zeros(m,2) ones(m,1) dUy*PP];
        Tz = [];
    end
    if(d == 3)
        Tx = [zeros(m,1) ones(m,1) zeros(m,2) dUx*PP];
        Ty = [zeros(m,2) ones(m,1) zeros(m,1) dUy*PP];
        Tz = [zeros(m,3) ones(m,1) dUz*PP];
    end
    
    kernel = PP'*K*PP;
    %get deformation structures
    optStruct.basis = basis;
    optStruct.kernel = kernel;%this is the factored kernel for tps, not for grbf!
    optStruct.K = K;
    optStruct.U = U;
    optStruct.Tx = Tx;
    optStruct.Ty = Ty;
    optStruct.Tz = Tz;
    x0 = [config.init_affine ...
        reshape(config.init_nonrigid', 1, numel(config.init_nonrigid))];
end

switch lower(config.normspec)
    %these are all bregmann now
    case 'restricted'
        
        init_flips = config.init_flips;
        init_tps = config.init_nonrigid;  % it should always be of size d*(n-d-1)
        
        if isempty(config.init_affine)
            % for your convenience, [] implies default affine
            config.init_affine = repmat([zeros(1,d) 1],1,d);
        end
        
        if config.opt_affine % optimize both affine and tps
            x0 = [x0 init_flips ];%init_pa];%OFF
        else % optimize tps only
            x0 = init_tps(end+1-d*(n-d-1):end);
            x0 = [x0 init_flips ];%init_pa];%OFF
        end
        
    case 'free'
        
        %init_pa = config.init_pa;%OFF
        if isempty(config.init_affine)
            % for your convenience, [] implies default affine
            config.init_affine = repmat([zeros(1,d) 1],1,d);
        end
        
        if config.opt_affine % optimize both affine and tps
            init_affine = [ ];
            x0 = [x0 0*config.Vtarget(:)'];
        else % optimize tps and normals
            
            x0 = [x0  0*config.Vtarget(:)'];
            
        end
        
    case 'fixed'
        %already done with this
end

optStruct.x0 = x0;

end