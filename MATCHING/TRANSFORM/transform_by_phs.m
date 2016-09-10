
function [warped_pts, bending_energy, vnew] = transform_by_phs(param, landmarks, ctrl_pts, ... 
   varargin)% Vmodel, kappa)
%this should basically be the same as by_tps since it only depends on the

if(length(varargin) == 2)
    Vmodel = varargin{1};
    kappa = varargin{2};
elseif(length(varargin) == 1)
    kappa = varargin{1};
end

%kernel
if (nargin==3)
    [n,d] = size(ctrl_pts);
    [m d] = size(landmarks);
    
    if(~exist('kappa', 'var'))
        kappa = d;
    end
    [K,U, dUx, dUy, dUz] = compute_kernel2(ctrl_pts, landmarks, 0, 'polyharmonic', kappa);
    %Pm = [ones(m,1) landmarks];
    %Pn = [ones(n,1) ctrl_pts];
    %PP = null(Pn'); 
    %B = [Pm U*PP]; 
    p2 = reshape(param(1:d*n),n,d);
    warped_pts = U*p2;
    tps_param =  reshape(param(d*(d+1)+1:d*n),n-d-1, d);
    bending_energy = trace(tps_param'*tps_param);
    vnew = [];
else
    [m,d] = size(landmarks);
    [n,d] = size(ctrl_pts);
    if(~exist('kappa', 'var'))
        kappa = d;
    end
    [K,U, dUx, dUy, dUz] = compute_kernel2(ctrl_pts,landmarks, 0, 'polyharmonic', kappa);
    %Pm = [ones(m,1) landmarks];
    %Pn = [ones(n,1) ctrl_pts];
    %PP = null(Pn'); 
    %B = [Pm U*PP]; 
    
    
    p2 = reshape(param(1:d*n),d,n)';
    warped_pts = U*p2;

        
        if(exist('Vmodel', 'var'))
            if(d == 2)
                Tx = dUx;%[zeros(m,1) ones(m,1) zeros(m,1) dUx*PP];
                Ty = dUy;%[zeros(m,2) ones(m,1) dUy*PP];
                T1 = Tx*p2;
                T2 = Ty*p2;
            end
            if( d == 3)
                Tx = [zeros(m,1) ones(m,1) zeros(m,2) dUx*PP];
                Ty = [zeros(m,2) ones(m,1) zeros(m,1) dUy*PP];
                Tz = [zeros(m,3) ones(m,1) dUz*PP];
                T1 = Tx*p2;
                T2 = Ty*p2;
                T3 = Tz*p2;
            end
            
        if( d == 2)
            for i = 1:size(T1,1)
                    Ai = [T1(i,:); T2(i,:)]';



                %[U1 S1 V1] = svd(Ai);
                %Rot = (U1*V1');
                Vd(i,:) = (Ai*Vmodel(d*(i-1)+1:d*i));
            end
            Vlong = reshape(Vd', 1, d*m)';
        %{
            %}
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
        
    tps_param = reshape(param(d*(d+1)+1:d*n)', n-d-1, d);
    bending_energy = trace(tps_param'*PP'*K*PP*tps_param);
end

