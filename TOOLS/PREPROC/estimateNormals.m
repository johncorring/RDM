function [v1 v2 v3] = estimateNormals(x, k, ordered)

    if(~exist('ordered', 'var' ))
        ordered = 0;
    end

    [n d] = size(x);

    v1 = [];
    v2 = [];
    v3 = [];

    if (d == 2)
        v1 = estimate2dNormals(x, k, ordered);
    elseif (d == 3)
        %%figure out how to input mesh ordering
        [theta phi] = est3dNormals(x, k);
        [v1 v2 v3] = sph2cart(theta, phi, ones(length(theta),1));
    else
        error('First argument should be an n x d matrix, with d = 2 or 3'); 
    end

end
