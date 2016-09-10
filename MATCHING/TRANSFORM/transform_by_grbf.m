
function [warped_pts, bending_energy, vnew] = transform_by_grbf(param, landmarks, ctrl_pts, Vmodel, flag)

[n,d] = size(ctrl_pts);


grbf_param = reshape(param(1:d*n),d,n);

[m,d] = size(landmarks);
[n,d] = size(ctrl_pts);
[K,U, dUx, dUy, dUz] = compute_kernel2(ctrl_pts,landmarks, 0, 'grbf');

grbf_param = reshape(param(1:d*n),d,n);

warped_pts = landmarks + U*grbf_param';
if(exist('Vmodel', 'var'))
    
    bending_energy = 0;%sum(sum(grbf_param'*K*grbf_param));
    
    
    if(d == 2)
        Tx = [dUx];
        Ty = [dUy];
        T1 = Tx*grbf_param';
        T2 = Ty*grbf_param';
    end
    if( d == 3)
        Tx = [dUx];
        Ty = [dUy];
        Tz = [dUz];
        T1 = Tx*grbf_param';
        T2 = Ty*grbf_param';
        T3 = Tz*grbf_param';
    end
    
    if( d == 2)
        for i = 1:size(T1,1)
            Ai = eye(2) + [T1(i,:); T2(i,:)]';
            
            
            
            %[U1 S1 V1] = svd(Ai);
            %Rot = (U1*V1');
            Vd(i,:) = (Ai*Vmodel(d*(i-1)+1:d*i));
        end
        Vlong = reshape(Vd', 1, d*m)';
        %{
            AA(:,1,:) = T1;
            AA(:,2,:) = T2;

            vec1 = AA(:,:,1)';
            vec2 = AA(:,:,2)';

            Vd(:,1) = vec1(:).*Vmodel;
            Vd(:,2) = vec2(:).*Vmodel;
            Vd = sum(Vd,2);
            Vlong = reshape(Vd', 1, d*m)';
            Vd = reshape(Vlong', d, m)';
        %}
        %vnew = Vtarget.*kron(flips_param',[1; 1; 1]);
        %Vlong = reshape(Vd', 1, d*m)';
        
    elseif( d == 3)
        for i = 1:size(T1,1)
            
            Ai = [T1(i,:); T2(i,:); T3(i,:)]';
            Vd(i,:) = (Ai*Vmodel(d*(i-1)+1:d*i));
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
        Vlong = reshape(Vd', 1, d*m)';
        Vd = reshape(Vlong', d, m)';
    end
    vnew = Vlong;
    
    
end

