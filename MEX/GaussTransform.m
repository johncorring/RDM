%%  $Author: bjian $
function [f,g] = GaussTransform(A,B,scale)	

if exist('mex_GaussTransform','file')
    [f,g] = mex_GaussTransform(A',B',scale);
    g = g';
else
    message = ['Precompiled GaussTransform module not found.\n' ...
        'If the corresponding MEX-functions exist, run the following command:\n' ...
        'mex mex_GaussTransform.c GaussTransform.c -output mex_GaussTransform'];
    message_id = 'MATLAB:MEXNotFound';
    error (message_id, message);
    
end
