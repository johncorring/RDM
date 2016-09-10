function [energy, grad] = rdmL2Bregmann_Free(param, optStruct)


    [temp_after_nonrigid, normals_after_nonrigid, bending, nonRigidParam] = ...
        localTransform(param, optStruct);
    Vd = reshape(normals_after_nonrigid', optStruct.d, ... 
        numel(normals_after_nonrigid)/optStruct.d)';
    VTargetUpdate = optStruct.Vtarget + param(end - length(optStruct.Vtarget) ...
        +1:end)';
    
    [energy, gradm2, gradv2, gradpa2] = ...
    gab_corr_costfunct(temp_after_nonrigid, optStruct.target, ... 
    normals_after_nonrigid, VTargetUpdate, optStruct.tempWts, ...
    optStruct.targWts, optStruct.scale, optStruct.lambda);

    [energy, gradm1, gradv1, gradpa1] = ...
    gab_corr_costfunct(optStruct.target, temp_after_nonrigid, ... 
    VTargetUpdate, normals_after_nonrigid, optStruct.targWts, ...
    optStruct.tempWts, optStruct.scale, optStruct.lambda);
    
    %compute energy plus regularization term
    energy = energy + optStruct.beta * bending;
    %gradient of tps
    grad =  (optStruct.basis'*gradm2)  ;

    %contributions to gradient of transform from the jacobians
    Tx = optStruct.Tx;
    Ty = optStruct.Ty;
    Tz = optStruct.Tz;
    d = optStruct.d;
    for j = 1:d
        grad(:,j) = grad(:,j) + (sum(Tx.*repmat(Vd(:,j).*gradv2(:,1),1,size(Tx,2)),1)');
        grad(:,j) = grad(:,j) + (sum(Ty.*repmat(Vd(:,j).*gradv2(:,2),1,size(Ty,2)),1)');
        if(d == 3)
            grad(:,j) = grad(:,j) + (sum(Tz.*repmat(Vd(:,j).*gradv2(:,3),1,size(Tz,2)),1)');
        end
    end
    
    if(strcmp(optStruct.motion, 'tps'))
        grad(d+2:end,:) = grad(d+2:end,:) + 2*optStruct.beta*optStruct.kernel*(nonRigidParam - optStruct.oldP')';
    elseif(strcmp(optStruct.motion, 'grbf'))
        grad = grad + 2*optStruct.beta*optStruct.kernel*(nonRigidParam - optStruct.oldP')';
        
    end

    grad = [grad; gradv1];

    grad = reshape(grad',1,numel(grad));
end


function [f g gv gpb] = gab_corr_costfunct(A, B, Va, Vb, Pa, Pb, scale, lambda)
    [f1 g1 gv1 gpb1] = GaborTransform(A, A, Va, Va, Pa, Pa, lambda, scale);
    [f2 g2 gv2 gpb2] = GaborTransform(A, B, Va, Vb, Pa, Pb, lambda, scale);
    [f3] = GaborTransform(B, B, Vb, Vb, Pb, Pb, lambda, scale);
    f = 2-2*f2/sqrt((f1)*(f3));
    g = -(g2 * sqrt(f1 * f3) - g1 * f2*f3/sqrt(f1*f3))/(f1* f3);%2*g1 - 2*g2;
    gv = -(gv2 * sqrt(f1 * f3) - gv1 * f2*f3/sqrt(f1*f3))/(f1* f3);
    gpb = -(gpb2* sqrt(f1 * f3) - gpb1 * f2*f3)/(f1^2 * f3^2);
end

function[temp_after_nonrigid, normals_after_nonrigid, bending, nonRigidParam] = ...
    localTransform(param, optStruct)
d = optStruct.d;
n = optStruct.n;
oldP = optStruct.oldP;

if(strcmp(optStruct.motion , 'tps'))

    affine_param = reshape(param(1:d*(d+1)),d,d+1);
    affine_param = affine_param';
    tps_param = reshape(param(d*(d+1)+1:d*n),d,n-d-1);
	nonRigidParam = tps_param;
    %apply tps. NB: QR is used to generate basis
    temp_after_nonrigid = optStruct.basis*[affine_param;tps_param'];
    bending = trace((tps_param - oldP')*optStruct.kernel*(tps_param - oldP')');

    %compute jacobians. TODO: untangle these and make the loop faster
    T1 = optStruct.Tx*[affine_param;tps_param'];
    T2 = optStruct.Ty*[affine_param;tps_param'];
    if(d == 3)
        T3 = optStruct.Tz*[affine_param; tps_param'];
    end

    %size of the template
    m = size(T1,1);



    if( d == 2)
        for i = 1:size(T1,1)
            Ai = [T1(i,:); T2(i,:)]';
            Vd(i,:) = (Ai*optStruct.Vtemplate(d*(i-1)+1:d*i));
        end

        %{
        AA(:,1,:) = T1;
        AA(:,2,:) = T2;

        vec1 = (AA(:,:,1)');
        vec2 = (AA(:,:,2)');

        Vd(:,1) = vec1(:).*Vtemplate;
        Vd(:,2) = vec2(:).*Vtemplate;
        Vd = sum(Vd,2);
        Vlong = reshape(Vd', 1, d*m)';
        Vd = reshape(Vlong', d, m)';
        %}
        %vnew = Vtarget.*kron(flips_param',[1; 1; 1]);
        %Vlong = reshape(Vd', 1, d*m)';
        normals_after_nonrigid = reshape(Vd', 1, d*m)';  

    elseif( d == 3)
        for i = 1:size(T1,1)

            Ai = [T1(i,:); T2(i,:); T3(i,:)]';
            Vd(i,:) = (Ai*optStruct.Vtemplate(d*(i-1)+1:d*i));
        end
            %{
        AA(:,1,:) = T1;
        AA(:,2,:) = T2;
        AA(:,3,:) = T3;

        Vd(:,1) = vec(AA(:,:,1)').*Vtemplate;
        Vd(:,2) = vec(AA(:,:,2)').*Vtemplate;
        Vd(:,3) = vec(AA(:,:,3)').*Vtemplate;
        Vlong = sum(Vd,2);
        %}
        normals_after_nonrigid = reshape(Vd', 1, d*m)';     
        %Vd = reshape(Vlong', d, m)';
    end

    %vnew = Vtarget;% + param(end-length(Vtarget)+1:end)';

elseif(strcmp(optStruct.motion , 'grbf'))

    grbf_param = reshape(param(1:d*size(optStruct.K,2)),d,size(optStruct.K,2));
    nonRigidParam = grbf_param;
    
    temp_after_nonrigid = optStruct.template + optStruct.U*grbf_param';
    bending = trace((grbf_param-oldP')*optStruct.K*(grbf_param-oldP')');

    %compute jacobians. TODO: untangle these and make the loop faster

    T1 = optStruct.Tx*grbf_param';
    T2 = optStruct.Ty*grbf_param';
    if(d == 3)
        T3 = optStruct.Tz*grbf_param';
    end

    %size of the template
    m = size(T1,1);
    Vd = zeros(size(T1));
    if( d == 2)
        for i = 1:size(T1,1)
                Ai = eye(2) + [T1(i,:); T2(i,:)]';
                Vd(i,:) = (Ai*optStruct.Vtemplate(d*(i-1)+1:d*i));
        end
        normals_after_nonrigid = reshape(Vd', 1, d*m)';
    elseif( d == 3)
        for i = 1:size(T1,1)
            Ai = [T1(i,:); T2(i,:); T3(i,:)]';
            Vd(i,:) = (Ai*optStruct.Vtemplate(d*(i-1)+1:d*i));
        end
        normals_after_nonrigid = reshape(Vd', 1, d*m)';     
    end
    %vnew = Vtarget + param(end-length(Vtarget)+1:end)';
end
   
end