function [energy, grad] = rdmL2Free_Matrix(param, init_affine, basis, kernel, ... 
    Tx, Ty, Tz, Vtarget, target, Vtemplate, scale, lambda, alpha, beta, n, d)
%Written by John Corring; email me at johncorring@gmail.com. RDM 2015 v1.


%get number of points in target
[na] = length(Vtarget);
n1 = na/d;

%where does the tps stop?
stN = d*size(Tx,2);

%figure out affine optimization strategy
if isempty(init_affine)
    %coeff_param = param(end-na/d+1:end);
    affine_param = reshape(param(1:d*(d+1)),d,d+1);
    affine_param = affine_param';
    tps_param = reshape(param(d*(d+1)+1:d*n),d,n-d-1);
    tps_param = tps_param';
else
    
    %coeff_param = param(end-na/d+1:end);
    stN = d*(n-d-1);
    tps_param = reshape(param(1:stN),d,n-d-1);
    tps_param = tps_param';
    affine_param = reshape(init_affine,d,d+1);
    affine_param = affine_param';
end

%apply tps. NB: QR is used to generate basis
after_tps = basis*[affine_param;tps_param];
bending = trace(tps_param'*kernel*tps_param);

%compute jacobians. TODO: untangle these and make the loop faster
T1 = Tx*[affine_param;tps_param];
T2 = Ty*[affine_param;tps_param];
if(d == 3)
    T3 = Tz*[affine_param; tps_param];
end

%size of the template
m = size(T1,1);

MatX = zeros(2*size(T1,1), 2);
MatY = zeros(2*n1, 2);

for b = 1:n1
    MatY(d*(b-1)+1:d*b,:) = eye(d)/scale^2;

end
if( d == 2)
    for i = 1:size(T1,1)
            Ai = [T1(i,:); T2(i,:)]';

            %[U1 S1 V1] = svd(Ai);
            %Rot = (U1*V1');
            Vd(i,:) = (Ai*Vtemplate(d*(i-1)+1:d*i));
            
            Rot = [Vd(i,:)/norm(Vd(i,:)); [-Vd(i,2) Vd(i,1)]/norm(Vd(i,:))];
           
            MatX(2*(i-1)+1:2*(i-1)+2,:) = Rot' * [1.2/scale^2 0 ; ...
                0 .8/scale^2] * Rot;
            %MatX(2*(i-1)+1:2*(i-1)+2,:) = MatX(2*(i-1)+1:2*(i-1)+2,:)/...
            %    norm(MatX(2*(i-1)+1:2*(i-1)+2,:), 'fro');
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

elseif( d == 3)
    for i = 1:size(T1,1)

        Ai = [T1(i,:); T2(i,:); T3(i,:)]';
        Vd(i,:) = (Ai*Vtemplate(d*(i-1)+1:d*i));
        
        [Rot b] = qr([Vd(i,:); eye(d)]');

        MatX(2*(i-1)+1:2*(i-1)+2,:) = Rot' * [1.2/scale^2 0 ; ...
            0 .8/scale^2] * Rot;
        %MatX(3*(i-1)+1:3*i,:) = eye(3)/scale^2;
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
    Vlong = reshape(Vd', 1, d*m)';     Vd = reshape(Vlong', d, m)';
end
       
vnew = Vtarget + param(end-length(Vtarget)+1:end)';

%reshaped template. unused
%Vq = reshape(Vtemplate',d,m)';

Vlong = reshape(Vd', 1, d*m)';

%compute the energy and gradients
%first run gets gradients for the new normals
[energy, gradm1, gradv1, gradpa1] = ...
    gab_corr_costfunct(target, after_tps, vnew, Vlong, MatY, MatX, ...
    ones(n1,1)/n1, ones(m,1)/m,  lambda);

%second run get gradients for the tps params
[energy, gradm2, gradv2, gradpa2] = ...
    gab_corr_costfunct(after_tps, target, Vlong, vnew, MatX, MatY, ...
    ones(m,1)/m, ones(n1,1)/n1, lambda);

%point rejection. unused
%kappa = 1/m;

%compute energy plus regularization term
energy = alpha * energy + beta * bending;% + kappa1*sum((coeff_param - 1/(n1/d)).^2);

%gradient of tps
grad = alpha*( (basis'*gradm2) );

%Rx = T1./(arrayfun((@(a,b,c,d) a*d-b*c), [T1 T2]));
%Ry = T2./(arrayfun((@(a,b,c,d) a*d-b*c), [T1 T2]));

%contributions to gradient of tps from the jacobians

for j = 1:d
    grad(:,j) = grad(:,j) + alpha*(sum(Tx.*repmat(Vd(:,j).*gradv2(:,1),1,n),1)');
    grad(:,j) = grad(:,j) + alpha*(sum(Ty.*repmat(Vd(:,j).*gradv2(:,2),1,n),1)');
    if(d == 3)
        grad(:,j) = grad(:,j) + alpha*(sum(Tz.*repmat(Vd(:,j).*gradv2(:,3),1,n),1)');
    end
end
%{
for j = 1:size(grad,1)
    for k = 1:size(grad,2)
        Ejk = zeros(size(grad));
        Ejk(j,k) = 1;
        grad(j,k) = grad(j,k) + alpha* sum(sum((Tx*Ejk)'*(repmat(Vd(:,1),1,2).*gradv2)));
        grad(j,k) = grad(j,k) + alpha* sum(sum((Ty*Ejk)'*(repmat(Vd(:,2),1,2).*gradv2)));
    end
end
%}
grad(d+2:n,:) = grad(d+2:n,:) + 2*beta*kernel*tps_param;
grad = [grad; zeros(size(gradv1))];

%testing 
gamma = 1;
grad(n+1:end,:) = gamma*gradv1;

if isempty(init_affine)
    grad = grad';
    grad = reshape(grad,1,d*n+d*n1);
    
    %if using point rejection need more space here. not usable
    %grad(stN+n1/d+1:end) = 0;gradpa1 + 2*kappa1 * (coeff_param' - 1/(n1/d));%%OFF
else
    grad(1:d+1,:) = [ ];
    grad = grad';
    grad = reshape(grad,1,d*(n-d-1));
    
    %if using point rejection need more space here. not usable
    %grad(stN+n1/d+1:end) = 0;gradpa1;%%OFF
end

function [f g gv gpb] = gab_corr_costfunct(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda)
    [f1 g1 gv1 gpb1] = GaborTransform3(A, A, Va, Va, Ma, Ma, Pa, Pa, lambda);
    [f2 g2 gv2 gpb2] = GaborTransform3(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda);
    [f3] = GaborTransform3(B, B, Vb, Vb, Mb, Mb, Pb, Pb, lambda);
    f = 2-2*f2/sqrt((f1)*(f3));
    g = -(g2 * sqrt(f1 * f3) - g1 * f2*f3/sqrt(f1*f3))/(f1* f3);%2*g1 - 2*g2;
    gv = -(gv2 * sqrt(f1 * f3) - gv1 * f2*f3/sqrt(f1*f3))/(f1* f3);
    gpb = -(gpb2* sqrt(f1 * f3) - gpb1 * f2*f3)/(f1^2 * f3^2);

function [f g gv gpb] = gab_costfunct_mat(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda)
[f1 g1 gv1 gpb1] = GaborTransform3(A, A, Va, Va, Ma, Ma, Pa, Pa, lambda);
[f2 g2 gv2 gpb2] = GaborTransform3(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda);
[f3] = GaborTransform3(B, B, Vb, Vb, Mb, Mb, Pb, Pb, lambda);
f = f1 - 2*f2 + f3;
g = 2*g1 - 2*g2;
gv = 2*gv1 - 2*gv2;
gpb = 2*gpb1 - 2*gpb2;