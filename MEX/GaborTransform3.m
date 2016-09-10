function [f, g, gv, gpb] = GaborTransform3(A,B,Va, Vb, Ma, Mb, Pa, Pb, lambda)	
%the one with the matrices
if exist('mex_GaborTransform','file')
    %assume that the Ma, Mb consist of the inverse covariance matrices
    %(precision matrices) from the multivariate Gaussian perspective
    [f, g, gv, gpb] = mex_GaborTransform3(A',B', Va', Vb', Ma', Mb', Pa', Pb',  lambda);
    g = g';
    gv = gv';
    gpb = gpb'; %TODO: implement membership gradients
else
    message = ['Precompiled GaborTransform module not found.\n' ...
        'If the corresponding MEX-functions exist, run the following command:\n' ...
        'mex mex_GaborTransform.c GaborTransform.c -output mex_GaborTransform'];
    message_id = 'MATLAB:MEXNotFound';
    error (message_id, message);
    
end