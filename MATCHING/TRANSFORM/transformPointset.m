function [transformed_pointset, v2] = ... 
    transformPointset(parameter, config)

v2 = [];
motion = config.motion;

switch lower(motion)
    case 'tps'
        ctrl_pts = config.ctrl_pts;
        [transformed_pointset ee v2] = transform_by_tps(parameter, config.template, ... 
                ctrl_pts, config.Vtemplate, 1);
    case 'grbf'
        ctrl_pts = config.ctrl_pts;
        [transformed_pointset ee v2] = transform_by_grbf(parameter, config.template, ... 
                ctrl_pts, config.Vtemplate, 1);
    otherwise
        error('Unknown motion type');
end;


