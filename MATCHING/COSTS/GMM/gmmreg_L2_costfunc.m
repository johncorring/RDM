function [f,g] = gmmreg_L2_costfunc(param, config)
%%=====================================================================
%% $RCSfile: gmmreg_L2_costfunc.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================
model = config.template;
scene = config.target;
motion = config.motion;
scale = config.scale;
[transformed_model, Va] = transform_pointset(model, motion, param, config.Vtarget);
switch lower(config.motion)
    case 'rigid2d'
        [f, grad] = rigid_costfunc(transformed_model, scene, scale);
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        grad = grad*model;
        theta = param(3);
        r = [-sin(theta) -cos(theta);
             cos(theta)  -sin(theta)];
        g(3) = sum(sum(grad.*r));
    case 'rigid3d'
        
        [r,gq] = quaternion2rotation(param(1:4));
        grad = grad';
        gm = grad*model; 
        g(1) = sum(sum(gm.*gq{1}));
        g(2) = sum(sum(gm.*gq{2}));
        g(3) = sum(sum(gm.*gq{3}));
        g(4) = sum(sum(gm.*gq{4}));         
        g(5) = sum(grad(1,:));
        g(6) = sum(grad(2,:));
        g(7) = sum(grad(3,:));
    case 'affine2d'
        [f grad] = general_costfunc(transformed_model, scene, scale);
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        g(3:6) = reshape(grad(1:2,:)*model,1,4);
        
        
    case 'area-rigid'
        [f ] = area_costfunc(transformed_model, scene, scale);
    case 'area-affine'
        [f ] = area_costfunc(transformed_model, scene, scale);
        grad = zeros(size(transformed_model));
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        g(3:6) = reshape(grad(1:2,:)*model,1,4);
    case 'rigid2dwv'
        [n1 dim] = size(model);
        %[n2 dim] = size(config.scene);
        %[transformed_template, Vb] = transform_pointset(config.scene, motion, [param(1:3) ones(1,n2)], config.Vscene);
        %Vb = config.Vscene;
        %[transformed_model, Va] = transform_pointset(model, motion, param, config.Vmodel);
        %[f,grad,gradvNULL] = gab_costfunct(transformed_template, transformed_model, Vb, Va, config.scale, config.lambda);
        [f,gradm,gradv] = gab_costfunct(transformed_model, config.scene, Va, config.Vscene, config.scale, config.lambda);
        
        g = zeros(3+n1,1);
        gradm = gradm';
        gradv = gradv';
        gradvt = gradv(1:2,:)*reshape(config.Vmodel',2,n1)';
        
        g(1) =sum(gradm(1,:));
        g(2) =sum(gradm(2,:));
        
        gradm = gradm(1:2,:)*model;%TODO: Transform V derivatives
        theta = param(3);
        r = [-sin(theta) -cos(theta);
             cos(theta)  -sin(theta)];
        g(3) = sum(sum(gradm.*r + gradvt.*r));
        
        gradv_all = (sum(gradv.*reshape(config.Vmodel',2,n1),1));
        g(4:end) = gradv_all;
        
    case 'affine2dwv'
        [n1 dim] = size(model);
        v1p = mod(config.Vmodel + param(3),2*pi);%estimateNormals(x1p, 2, 1);%
        v1pk = kron(v1p, ones(2,1));
        v1pk = v1pk - kron(ones(n1,1), pi/2*[0; 1]);
        Va = cos(v1pk);
        
        [n2 dim] = size(config.scene);
        v2p = mod(config.Vscene,2*pi);
        v2pk = kron(v2p,ones(2,1));
        v2pk = v2pk - kron(ones(n2,1), pi/2*[0; 1]);
        Vb = cos(v2pk);
        
        [f, grad] = gab_costfunct(transformed_model, scene, Va, Vb, config.scale, config.lambda);%TODO: GRAD
        grad = grad';
        g = zeros(1,6);
        %g(3:6) = reshape(grad*model,1,4);
        
        %grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        g(3:6) = reshape(grad*model,1,4);
    case 'affine3d'
        [f,grad] = general_costfunc(transformed_model, scene, scale);
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        g(3) = sum(grad(3,:));
        g(4:12) = reshape(grad*model,1,9);
    otherwise
        error('Unknown motion type');
end;


function [f, g] = rigid_costfunc(A, B, scale)
[f, g] =  GaussTransform(A,B,scale);
f = -f; g = -g;

function [f g gv] = gab_costfunct(A, B, Va, Vb, scale, lambda)
[f1 g1 gv1] = GaborTransform(A, A, Va, Va, lambda, scale);
[f2 g2 gv2] = GaborTransform(A, B, Va, Vb, lambda, scale);
f = f1 - 2*f2;
g = 2*g1 - 2*g2;
gv = 2*gv1 - 2*gv2;

function [f, g] = general_costfunc(A, B, scale)
[f1, g1] = GaussTransform(A,A,scale);
[f2, g2] = GaussTransform(A,B,scale);
f =  f1 - 2*f2;
g = 2*g1 - 2*g2;

function [f] = area_costfunc(A,B,scale)

[f1] = AreaTransform(A,B,scale);%areaGaussv1(A, B, scale);%gAreaFunc2(A,B,scale);%
[f2, g2] = GaussTransform(A,B,scale);
f = f1 - .1* f2;
%g = - g2;