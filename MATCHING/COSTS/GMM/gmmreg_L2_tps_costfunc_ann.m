function [energy, grad] = gmmreg_L2_tps_costfunc_ann(param, init_affine, basis, kernel, Tx, Ty, Tz, Vmodel, model, Vscene, scale, lambda, alpha, beta, n, d)

kappa1 = 1/scale;
[na] = length(Vmodel);
n1 = na/d;
stN = d*size(Tx,2);
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
after_tps = basis*[affine_param;tps_param];
bending = trace(tps_param'*kernel*tps_param);

T1 = Tx*[affine_param;tps_param];
T2 = Ty*[affine_param;tps_param];
if(d == 3)
   T3 = Tz*[affine_param; tps_param]; 
end
m = size(T1,1);

for i = 1:size(T1,1)
if( d == 2)

    Ai = [T1(i,:); T2(i,:)]';

end
if( d == 3)

    Ai = [T1(i,:); T2(i,:); T3(3,:)]';
end
    [U1 S1 V1] = svd(Ai);
    Rot = (U1*V1');
    Vd(i,:) = (Ai*Vscene(d*(i-1)+1:d*i)); 
    %Vd(i,:) = Vscene(d*(i-1)+1:d*i);
end
vnew = Vmodel + param(d*n+1:end-1)';
Vq = reshape(Vscene',d,m)';

%}
%[Ug sg Vg] = svd(affine_param(2:3,:));
%Vd(:,1) = sum(T1.*Vq,2);
%Vd(:,2) = sum(T2.*Vq,2);
%Vd = Vd*(Ug*Vg')';
Vlong = reshape(Vd', 1, d*m)';

[energy, gradm1, gradv1, gradpa1] = ... 
    gab_costfunct(model, after_tps, vnew, Vlong, ones(n1,1)/n1, ones(m,1)/m,  scale, lambda);
[energy, gradm2, gradv2, gradpa2] = ... 
    gab_costfunct(after_tps, model, Vlong, vnew, ones(m,1)/m, ones(n1,1)/n1, scale, lambda);

energy = alpha * energy + beta * bending;% + kappa1*sum((coeff_param - 1/(n1/d)).^2);

grad = alpha*( (basis'*gradm2) );

%Rx = T1./(arrayfun((@(a,b,c,d) a*d-b*c), [T1 T2]));
%Ry = T2./(arrayfun((@(a,b,c,d) a*d-b*c), [T1 T2]));

for j = 1:d
grad(:,j) = grad(:,j) + alpha*(sum(Tx.*repmat(Vd(:,j).*gradv2(:,1),1,n),1)');
grad(:,j) = grad(:,j) + alpha*(sum(Ty.*repmat(Vd(:,j).*gradv2(:,2),1,n),1)');
end
%}
%grad(:,2) = grad(:,2) + alpha*(ones(n,m)*(Vd(:,2).*gradv2(:,2)));
if(d == 3)
    grad(:,3) = grad(:,3) + alpha*(Tz'*(sum(reshape(Vscene',d,m)',2).*gradv2(:,3)));
end


%grad = grad + alpha*([sum((Tx'*gradv2),2) sum(Ty'*gradv2,2)]);
grad(d+2:n,:) = grad(d+2:n,:) + 2*beta*kernel*tps_param;
grad = [grad; zeros(size(gradv1))];
grad(n+1:end,:) = alpha*gradv1;%reshape(gradv1', d, n1)';

if isempty(init_affine) 
    %% In this case, the length of gradient should be n*d    
    grad = grad';
    grad = reshape(grad,1,d*n+d*n1);
    grad = [grad];
    %grad = [grad zeros(1, 2*size(gradv1,1))];
    %grad(stN+1:stN+n1/d) = alpha*(sum(gradv1'.*reshape(Vmodel',d,n1/d),1));
    %grad(stN+n1/d+1:end) = 0;gradpa1 + 2*kappa1 * (coeff_param' - 1/(n1/d));%%OFF
else 
    %% In this case, the length of parameter should be (n-d-1)*d    
    grad(1:d+1,:) = [ ];
    grad = grad';
    grad = reshape(grad,1,d*(n-d-1));
    grad = [grad];
    %grad = [grad zeros(1, 2*size(gradv1,1))];
    %grad(stN+1:stN+n1/d) = (sum(gradv1'.*reshape(Vmodel',d,n1/d),1));
    %grad(stN+n1/d+1:end) = 0;gradpa1;%%OFF
end

function [f, g] = general_costfunc(A, B, scale)
[f1, g1] = GaussTransform(A,A,scale);
[f2, g2] = GaussTransform(A,B,scale);
f =  f1 - 2*f2;
g = 2*g1 - 2*g2;

function [f g gv gpb] = gab_costfunct(A, B, Va, Vb, Pa, Pb, scale, lambda)
[f1 g1 gv1 gpb1] = GaborTransform(A, A, Va, Va, Pa, Pa, lambda, scale);
[f2 g2 gv2 gpb2] = GaborTransform(A, B, Va, Vb, Pa, Pb, lambda, scale);
[f3] = GaborTransform(B, B, Vb, Vb, Pb, Pb, lambda, scale);
f = f1 - 2*f2 + f3;
g = 2*g1 - 2*g2;
gv = 2*gv1 - 2*gv2;
gpb = 2*gpb1 - 2*gpb2;