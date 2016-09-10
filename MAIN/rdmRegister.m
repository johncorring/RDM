
function [param, transformedTemplate, transformedNormal, history, config] = ... 
    rdmRegister(config)


history.x = [ ];
history.fval = [ ];

if nargin~=1
    error('Create a configuration object with rdmConfig and call as  rdmRegister(config).');
end

[n,d] = size(config.ctrl_pts); 

if (d~=2)&&(d~=3)
    error('This program only handles 2d or 3d pointsets in the format of [numsamples x numdimensions]');
end


options = optimset( 'display','iter', 'LargeScale','off','GradObj','on', ... 
    'TolFun',1e-15, 'TolX',1e-15, 'TolCon', 1e-15);

options = optimset(options, 'outputfcn', @outfun);
options = optimset(options, 'MaxIter', config.max_iter);
options = optimset(options, 'GradObj', 'on');

if(config.gradflag == 0)
    options = optimset(options, 'GradObj', 'off');
end

tic

optStruct = createOptStruct(config);
x0 = optStruct.x0;
[n,d] = size(config.ctrl_pts); 
[m1, d] = size(optStruct.target);

switch lower(config.normspec)
    case 'fixed'
        param = fminunc(@(x)rdmL2Bregmann_Fixed(x, optStruct), x0, options);
        [transformedTemplate, transformedNormal] = ...
            transformPointset(param, config);
        
        %unload parameters
        if(strcmp(config.motion, 'tps'))
                config.init_nonrigid = reshape(param(d*(d+1)+1:d*n), d, n-d-1)';
                config.init_affine = param(1:d*(d+1));;
            
        elseif(strcmp(config.motion, 'grbf'))
            config.init_nonrigid = reshape(param(1:d*n), d, n)';
        end
            config.oldP = config.init_nonrigid;
    case 'free'
        
        param = fminunc(@(x) rdmL2Bregmann_Free(x, optStruct), x0, options);
        [transformedTemplate, transformedNormal] = ...
            transformPointset(param, config);
        
        %unload parameters
        if config.opt_affine
            %config.init_pa = param(end-m1+1:end);
            %config.init_flips = param(end-2*m1+1:end-m1);
            config.Vtarget = param(end-d*m1+1:end)';
            if(strcmp(config.motion, 'tps'))
                config.init_nonrigid = reshape(param(d*(d+1)+1:d*n), d, n-d-1)';
                 config.init_affine = param(1:d*(d+1));
            
            elseif(strcmp(config.motion, 'grbf'))
                config.init_nonrigid = reshape(param(1:d*n), d, n)';
            end
            config.oldP = config.init_nonrigid;
            
        else
            config.init_nonrigid = param(1:end-2*m1);
            config.oldP = config.init_nonrigid;
            %config.init_flips = param(end-2*m1+1:end-m1);
            %config.init_pa = param(end-m1+1:end);
        end

    case 'restricted'
        
        param = fminunc(@(x) rdmL2Bregmann_Restricted(x, optStruct), x0, options);
        [transformedTemplate transformedNormal] = ...
            transformPointset(param, config);
        
        %unload parameters
        
        if config.opt_affine
            %config.init_pa = param(end-m1+1:end);
            %config.init_flips = param(end-2*m1+1:end-m1);
            config.Vtarget = param(end-d*m1+1:end)';
            config.init_nonrigid = param(end+1-d*(n-d-1)-2*m1:end-2*m1);
            config.init_affine = param(1:d*(d+1));
            config.oldP = config.init_nonrigid;
        else
            config.init_nonrigid = param(1:end-d*m1);
            config.init_flips = param(end-d*m1+1:end-m1);
            config.oldP = config.init_nonrigid;
            %config.init_pa = param(end-m1+1:end);
        end
        
    otherwise
        x0 = config.init_param;
        param = fmincon(@gmmreg_L2_costfunc, x0, [ ],[ ],[ ],[ ], config.Lb, config.Ub, [ ], options, config);
        [transformed_model v] = transform_pointset(config.template, config.motion, param, config.Vtemplate);
        config.init_param = param;
end
toc

function stop = outfun(x,optimValues,state)
    stop = false;

    switch state
        case 'init'
            hold on
        case 'iter'
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x];
        case 'done'
            hold off
        otherwise
    end
end

end
