% Perform a spatial tranformation on a given pointset
% motion:  the motion model represented by string, can be
%         'rigid2d',    'rigid3d', 'affine2d',  'affine3d', 'tps'
% parameter: a row vector
function [transformed_pointset, v2] = transform_pointset(pointset, motion, parameter, varargin)
%%=====================================================================
%% $RCSfile: transform_pointset.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================
v2 = [];
switch lower(motion)
    case 'rigid2d'
        transformed_pointset = transform_by_rigid2d(pointset, parameter);
    case 'rigid3d'
        transformed_pointset = transform_by_rigid3d(pointset, parameter);
    case 'area-affine'
        transformed_pointset = transform_by_affine2d(pointset, parameter);
    case 'affine2d'
        transformed_pointset = transform_by_affine2d(pointset, parameter);
    case 'rigid2dwv'
        v = varargin{1};
        [transformed_pointset v2] = transform_by_rigid2dwv(pointset, v, parameter);
    case 'affine2dwv'
        v = varargin{1};
        [transformed_pointset v2] = transform_by_affine2dwv(pointset, v,  parameter);
    case 'affine3d'
        transformed_pointset = transform_by_affine3d(pointset, parameter);
    case 'tps'
        ctrl_pts = varargin{1};
        init_affine = varargin{2};
        [n,d] = size(ctrl_pts);
        p = reshape([init_affine parameter],d,n); p = p'; 
        transformed_pointset = transform_by_tps_GMM(p, pointset, ctrl_pts);
    case 'tps-restricted'
        if(length(varargin) == 4)
                flag = varargin{4};
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
        else 

            if(length(varargin) == 3)
            ctrl_pts = varargin{1};
            init_affine = varargin{2};
            v = varargin{3};
            p = [init_affine(:)' parameter(:)'];
            flag = 1;
            else
                p = [parameter(:)];
                if(length(varargin) == 2)            
                ctrl_pts = varargin{1};
                v = varargin{2};

                else
                    if(length(varargin) == 1)
                        ctrl_pts = varargin{1};
                        flag = 0;
                    else
                        error('wrong number of varargin arguments');
                    end
                end

            end
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_tps(p, pointset, ctrl_pts, v, flag);
        else
            [transformed_pointset] = transform_by_tps(p, pointset, ctrl_pts);
        end
        
    case 'tps-free'
        if(length(varargin) == 4)
                flag = varargin{4};
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
        else 

            if(length(varargin) == 3)
            ctrl_pts = varargin{1};
            init_affine = varargin{2};
            v = varargin{3};
            p = [init_affine(:)' parameter(:)'];
            flag = 1;
            else
                p = [parameter(:)];
                if(length(varargin) == 2)            
                ctrl_pts = varargin{1};
                v = varargin{2};

                else
                    if(length(varargin) == 1)
                        ctrl_pts = varargin{1};
                        flag = 0;
                    else
                        error('wrong number of varargin arguments');
                    end
                end

            end
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_tps(p, pointset, ctrl_pts, v, flag);
        else
            [transformed_pointset] = transform_by_tps(p, pointset, ctrl_pts);
        end
          
    case 'tps-matrix'
        if(length(varargin) == 4)
                flag = varargin{4};
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
        else 

            if(length(varargin) == 3)
            ctrl_pts = varargin{1};
            init_affine = varargin{2};
            v = varargin{3};
            p = [init_affine(:)' parameter(:)'];
            flag = 1;
            else
                p = [parameter(:)];
                if(length(varargin) == 2)            
                ctrl_pts = varargin{1};
                v = varargin{2};

                else
                    if(length(varargin) == 1)
                        ctrl_pts = varargin{1};
                        flag = 0;
                    else
                        error('wrong number of varargin arguments');
                    end
                end

            end
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_tps(p, pointset, ctrl_pts, v, flag);
        else
            [transformed_pointset] = transform_by_tps(p, pointset, ctrl_pts);
        end
    case 'polyharmonic'
        if(length(varargin) == 4)
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
                kappa = varargin{4};
        elseif(length(varargin) == 3)
                %assume v was dropped
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                p = [init_affine(:)' parameter(:)'];
                kappa = varargin{3};
        elseif(length(varargin) == 2)
                %assume flag was dropped
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                p = [init_affine(:)' parameter(:)'];
                kappa = size(ctrl_pts,2);
        else
            error('wrong number of varargin arguments');
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_phs(p, pointset, ctrl_pts, v, kappa);
        else
            [transformed_pointset] = transform_by_phs(p, pointset, ctrl_pts, kappa);
        end
    case 'grbf-restricted'
        if(length(varargin) == 4)
                flag = varargin{4};
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
        else 

            if(length(varargin) == 3)
            ctrl_pts = varargin{1};
            init_affine = varargin{2};
            v = varargin{3};
            p = [init_affine(:)' parameter(:)'];
            flag = 1;
            else
                p = [parameter(:)];
                if(length(varargin) == 2)            
                ctrl_pts = varargin{1};
                v = varargin{2};

                else
                    if(length(varargin) == 1)
                        ctrl_pts = varargin{1};
                        flag = 0;
                    else
                        error('wrong number of varargin arguments');
                    end
                end

            end
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_grbf(p, pointset, ctrl_pts, v, flag);
        else
            [transformed_pointset] = transform_by_grbf(p, pointset, ctrl_pts);
        end
    case 'clamped'
        
        if(length(varargin) == 4)
                flag = varargin{4};
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
        else 

            if(length(varargin) == 3)
            ctrl_pts = varargin{1};
            init_affine = varargin{2};
            v = varargin{3};
            p = [init_affine(:)' parameter(:)'];
            flag = 1;
            else
                p = [parameter(:)];
                if(length(varargin) == 2)            
                ctrl_pts = varargin{1};
                v = varargin{2};

                else
                    if(length(varargin) == 1)
                        ctrl_pts = varargin{1};
                        flag = 0;
                    else
                        error('wrong number of varargin arguments');
                    end
                end

            end
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_clamped(p, pointset, ctrl_pts, v, flag);
        else
            [transformed_pointset] = transform_by_clamped(p, pointset, ctrl_pts);
        end
    case 'grbf-free'
        if(length(varargin) == 4)
                flag = varargin{4};
                ctrl_pts = varargin{1};
                init_affine = varargin{2};
                v = varargin{3};
                p = [init_affine(:)' parameter(:)'];
        else 

            if(length(varargin) == 3)
            ctrl_pts = varargin{1};
            init_affine = varargin{2};
            v = varargin{3};
            p = [init_affine(:)' parameter(:)'];
            flag = 1;
            else
                p = [parameter(:)];
                if(length(varargin) == 2)            
                ctrl_pts = varargin{1};
                v = varargin{2};

                else
                    if(length(varargin) == 1)
                        ctrl_pts = varargin{1};
                        flag = 0;
                    else
                        error('wrong number of varargin arguments');
                    end
                end

            end
        end
        %[n,d] = size(ctrl_pts);%reshape([init_affine parameter],d,n); p = p'; 
        %p = reshape(p,d,length(p)/d)';
        if(exist('v', 'var'))
            [transformed_pointset ee v2] = transform_by_grbf(p, pointset, ctrl_pts, v, flag);
        else
            [transformed_pointset] = transform_by_grbf(p, pointset, ctrl_pts);
        end
    otherwise
        error('Unknown motion type');
end;


