function [f, g, gv, gpb] = GaborTransform3_mat(A,B,Va, Vb, Ma, Mb, Pa, Pb, lambda)	
%the one with the matrices
    [f, g, gv, gpb] = GaborMat(A,B, Va, Vb, Ma, Mb, Pa, Pb, lambda);
    g = g;
    gv = gv;
    gpb = gpb; %TODO: implement gradient
end

function [f g gv gpb] = GaborMat(mA, mB, vA, vB, A, B, pA, pB, lambda)

[n1 d] = size(mA);
[n2 d] = size(mB);

g = zeros(n1, d);
gv = zeros(n1, d);
gpb = zeros(n1, 1);

f = 0;

%for each pair
for j = 1:n1
    %Amat is matrix for A set
    Amat = A((j-1)*d + 1:(j-1)*d+d,:);
    maj = mA(j,:);
    vaj = vA(j,:);
    
    
    for k = 1:n2
        %Bmat is matrix for B set
        Bmat = B((k-1)*d + 1 : (k-1)*d+d,:);
        mbk = mB(k,:);
        vbk = vB(k,:);
        
        %sum matrix
        sumMat = Amat + Bmat;
        
        qForm = (((vaj - vbk)/lambda - sqrt(-1) * (maj*Amat + mbk*Bmat)) * inv(sumMat)) * ...
            ((vaj - vbk)/lambda - sqrt(-1) * (maj*Amat + mbk*Bmat)).';
        
        dqFormdm = - 2.0 * sqrt(-1) *  (inv(sumMat) * ...
            ((vaj - vbk)/lambda - sqrt(-1) * (maj*Amat + mbk*Bmat)).').' * Amat;
        
        dqFormdv = 2.0 * (1.0 / lambda) * (inv(sumMat) * ...
            ((vaj - vbk)/lambda - sqrt(-1) * (maj*Amat + mbk*Bmat)).').';
        
        oPart = - maj * Amat * maj.' - mbk * Bmat * mbk.' - 2.0*sqrt(-1)/lambda * vaj * ...
            maj.' + 2.0*sqrt(-1)/lambda * vbk * mbk.';
        
        doPartdm = - 2.0 *  maj * Amat - 2.0 * sqrt(-1)/lambda * vaj;
        
        doPartdv = -2.0 * sqrt(-1)/lambda * maj;
        
        f = f + ((2*pi)^d/det(sumMat)).^.5 * pA(j) * pB(k) * exp( - .5 * qForm + .5 * oPart);
        g(j,:) = g(j,:) + (-.5*dqFormdm + .5 * doPartdm)* ... 
           ((2*pi)^d/det(sumMat)).^.5 * pA(j) * pB(k) * exp( - .5 * qForm + .5 * oPart);
        gv(j,:) = gv(j,:) + (-.5*dqFormdv + .5 * doPartdv)* ... 
           ((2*pi)^d/det(sumMat)).^.5 * pA(j) * pB(k) * exp( - .5 * qForm + .5 * oPart);
        
        
        gpb(j,:) = gpb(j,:) + ((2*pi)^d/det(sumMat)).^.5 * pB(k) * ... 
            exp( - .5 * qForm + .5 * oPart);
        
    end
    
end

f = real(f);
g = real(g);
gv = real(gv);
gpb = real(gpb);
end