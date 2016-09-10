function [f, g, gv, gpb] = GaborTransform(A,B,Va, Vb, Pa, Pb, lambda, scale)	
%the one with real arithmetic in C
if exist('mex_GaborTransform','file')
    [f, g, gv, gpb] = mex_GaborTransform(A',B', Va', Vb', Pa', Pb',  lambda, scale);
    g = g';
    gv = gv';
    gpb = gpb'; %TODO: implement gradient
else
    message = ['Precompiled GaborTransform module not found.\n' ...
        'If the corresponding MEX-functions exist, run the following command:\n' ...
        'mex mex_GaborTransform.c GaborTransform.c -output mex_GaborTransform'];
    message_id = 'MATLAB:MEXNotFound';
    error (message_id, message);
    
end