

function [energy, grad] = rdmL2Bregmann_Restricted(param, optStruct)
%TO IMPLEMENT
 

end

function rdmL2Restricted(param, init_affine, basis, kernel, Tx, Ty, Tz, Vtarget, target, Vtemplate, scale, lambda, alpha, beta, n, d)

[n1] = length(Vtarget);
stN = d*size(Tx,2);

if isempty(init_affine)
    flips_param = param(end-n1/d+1:end);
    %coeff_param = param(end-n1/d+1:end);
    affine_param = reshape(param(1:d*(d+1)),d,d+1);
    affine_param = affine_param';
    tps_param = reshape(param(d*(d+1)+1:d*n),d,n-d-1);
    tps_param = tps_param';
else
    flips_param = param(end-n1/d+1:end);
    %coeff_param = param(end-n1/d+1:end);
    stN = d*(n-d-1);
    tps_param = reshape(param(1:stN),d,n-d-1);
    tps_param = tps_param';
    affine_param = reshape(init_affine,d,d+1);
    affine_param = affine_param';
end

after_tps = basis*[affine_param;tps_param];
bending = trace(tps_param'*kernel*tps_param);

T1 = Tx*[affine_param;tps_param];
T2 = Ty*[affine_param;tps_param];
if(d == 3)
    T3 = Tz*[affine_param; tps_param];
end

m = size(T1,1);


if( d == 2)
    for i = 1:size(T1,1)
            Ai = [T1(i,:); T2(i,:)]';



        %[U1 S1 V1] = svd(Ai);
        %Rot = (U1*V1');
        Vd(i,:) = (Ai*Vtemplate(d*(i-1)+1:d*i));
    end
%{
    
     AA(:,1,:) = T1;
    AA(:,2,:) = T2;
    
    vec1 = AA(:,:,1)';
    vec2 = AA(:,:,2)';
    
    Vd(:,1) = vec1(:).*Vtemplate;
    Vd(:,2) = vec2(:).*Vtemplate;
    Vd = sum(Vd,2);
    Vlong = reshape(Vd', 1, d*m)';
    Vd = reshape(Vlong', d, m)';
    %}
    Vlong = reshape(Vd', 1, d*m)';
    %vnew = Vtarget.*kron(flips_param', [1; 1]);

    vnew = Vtarget.*kron(flips_param',[1; 1]);
elseif( d == 3)
        %{
    for i = 1:size(T1,1)

         Ai = [T1(i,:); T2(i,:); T3(i,:)]';
        Vd(i,:) = (Ai*Vtemplate(d*(i-1)+1:d*i));
    end
    %}
    AA(:,1,:) = T1;
    AA(:,2,:) = T2;
    AA(:,3,:) = T3;
    
    Vd(:,1) = vec(AA(:,:,1)').*Vtemplate;
    Vd(:,2) = vec(AA(:,:,2)').*Vtemplate;
    Vd(:,3) = vec(AA(:,:,3)').*Vtemplate;
    Vlong = sum(Vd,2);
    Vd = reshape(Vlong', d, m)';
    vnew = Vtarget.*kron(flips_param',[1; 1; 1]);
end
        

[energy, gradm1, gradv1, gradpa1] = ...
    gab_costfunct(target, after_tps, vnew, Vlong, ones(n1/d,1)/(n1/d), ones(m,1)/m,  scale, lambda);
[energy, gradm2, gradv2, gradpa2] = ...
    gab_costfunct(after_tps, target, Vlong, vnew, ones(m,1)/m, ones(n1/d,1)/(n1/d), scale, lambda);

energy = alpha*energy + beta * bending;% + kappa1*sum((coeff_param - 1/(n1/d)).^2);

grad = alpha*( (basis'*gradm2) );

for j = 1:d
    grad(:,j) = grad(:,j) + alpha*(sum(Tx.*repmat(Vd(:,j).*gradv2(:,1),1,n),1)');
    grad(:,j) = grad(:,j) + alpha*(sum(Ty.*repmat(Vd(:,j).*gradv2(:,2),1,n),1)');
    if(d == 3)
        grad(:,j) = grad(:,j) + alpha*(sum(Tz.*repmat(Vd(:,j).*gradv2(:,3),1,n),1)');
    end
end

grad(d+2:n,:) = grad(d+2:n,:) + beta*kernel*tps_param;


gamma = 0;
%grad(n+1:end,:) = gamma*gradv1;
if isempty(init_affine)
    %% In this case, the length of gradient should be n*d
    grad = grad';
    grad = reshape(grad,1,d*n);
    %grad = [grad zeros(1, size(gradv1,1))];
    grad = [grad alpha*gamma*(sum(gradv1'.*reshape(Vtarget',d,n1/d),1))];
    %grad(stN+n1/d+1:end) = 0;gradpa1 + 2*kappa1 * (coeff_param' - 1/(n1/d));%%OFF
else
    %% In this case, the length of parameter should be (n-d-1)*d
    grad(1:d+1,:) = [ ];
    grad = grad';
    grad = reshape(grad,1,d*(n-d-1));
    grad = [grad zeros(1, size(gradv1,1))];
    grad(stN+1:stN+n1/d) = (sum(gradv1'.*reshape(Vtarget',d,n1/d),1));
    %grad(stN+n1/d+1:end) = 0;gradpa1;%%OFF
end

function [f g gv gpb] = gab_costfunct(A, B, Va, Vb, Pa, Pb, scale, lambda)
[f1 g1 gv1 gpb1] = GaborTransform(A, A, Va, Va, Pa, Pa, lambda, scale);
[f2 g2 gv2 gpb2] = GaborTransform(A, B, Va, Vb, Pa, Pb, lambda, scale);
[f3] = GaborTransform(B, B, Vb, Vb, Pb, Pb, lambda, scale);
f = f1 - 2*f2 + f3;
g = 2*g1 - 2*g2;
gv = 2*gv1 - 2*gv2;
gpb = 2*gpb1 - 2*gpb2;