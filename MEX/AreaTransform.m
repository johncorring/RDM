function [f,g] = AreaTransform(A,B,scale)	

if exist('mex_GaussTransform','file')
    [f,g] = mex_AreaTransform(A',B',scale);
    g = g';
else
    message = ['Precompiled AreaTransform module not found.\n' ...
        'If the corresponding MEX-functions exist, run the following command:\n' ...
        'mex areaTransform.c -output areaTransform'];
    message_id = 'MATLAB:MEXNotFound';
    error (message_id, message);
    
end