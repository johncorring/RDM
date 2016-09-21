function [nonrigid affine] = getCorrsTransformation(cPts, ...
    X, Y, bendPenalty)
    
    [PP,kernel] = tps_set_ctrl_pts(cPts);
    [U,Pm,Q1,Q2,R] = tps_set_landmarks(X,cPts);
    params = tps_compute_param(PP,kernel,U,Pm,Q1,Q2,R,bendPenalty,Y)
    affine = params(1:4,1:3);
    nonrigid = params(5:end,:);
end
