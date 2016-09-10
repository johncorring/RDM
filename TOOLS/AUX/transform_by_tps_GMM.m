function [warped_pts, bending_energy] = transform_by_tps_GMM(param, landmarks, ctrl_pts)
%%=====================================================================
%% $Author: bing.jian $
%%=====================================================================
if (nargin==2)
    [n,d] = size(landmarks);
    [B,lambda] = compute_basis(landmarks);
    warped_pts = B*param;
    tps_param = param(d+2:n,:);
    bending_energy = trace(tps_param'*diag(lambda)*tps_param);
else
    [m,d] = size(landmarks);
    [n,d] = size(ctrl_pts);
    [K,U] = compute_kernel(ctrl_pts,landmarks);
    Pm = [ones(m,1) landmarks];
    Pn = [ones(n,1) ctrl_pts];
    PP = null(Pn'); B = [Pm U*PP]; 
    warped_pts = B*param;
    tps_param = param(d+2:n,:);
    bending_energy = trace(tps_param'*PP'*K*PP*tps_param);
end

