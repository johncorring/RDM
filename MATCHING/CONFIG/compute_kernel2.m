function [K,U, dUx, dUy, dUz] = compute_kernel2(ctrl_pts, landmarks, ...
    lambda, type, kappa)


[n,d] = size(ctrl_pts);

if(~exist('kappa', 'var'))
    kappa = d;
end

if (nargin<3)
    lambda = 0;
end

K = lambda*ones(n); %% K only depends on ctrl_pts and lambda

switch(type)
    case('tps')
        switch d
            case 2
                for i=1:n
                    for j=1:n
                        r = norm(ctrl_pts(i,1:2) - ctrl_pts(j,1:2));
                        if (r>0)
                            K(i,j) =   r*r*log(r);
                        end
                    end
                end
            case 3
                for i=1:n
                    for j=1:n
                        r = norm(ctrl_pts(i,1:3) - ctrl_pts(j,1:3));
                        K(i,j) =   -r;
                    end
                end
        end
        
        if (nargin>=2)
            [m,d] = size(landmarks);
            %U = zeros(m,n);
            U = lambda*ones(m,n);
            
            dUx = zeros(m,n);
            dUy = zeros(m,n);
            dUz = zeros(m,n);
            
            switch d
                case 2
                    for i=1:m
                        for j=1:n
                            r = norm(landmarks(i,1:2) - ctrl_pts(j,1:2));
                            if (r>0)
                                U(i,j) =  r^2*log(r);
                                dUx(i,j) = (landmarks(i,1)-ctrl_pts(j,1))*(2*log(r) +1);
                                dUy(i,j) = (landmarks(i,2)-ctrl_pts(j,2))*(2*log(r) +1);
                            end
                        end
                    end
                case 3
                    for i=1:m
                        for j=1:n
                            r = norm(landmarks(i,1:3) - ctrl_pts(j,1:3));
                            U(i,j) =   -r;
                            if(r > 0)
                                dUx(i,j) = - (landmarks(i,1) - ctrl_pts(j,1))/r;
                                dUy(i,j) = - (landmarks(i,2) - ctrl_pts(j,2))/r;
                                dUz(i,j) = - (landmarks(i,3) - ctrl_pts(j,3))/r;
                            end
                        end
                    end
            end
        end
    case('polyharmonic')
        %Work in progress
        for i=1:n
            for j=1:n
                r = norm(ctrl_pts(i,:) - ctrl_pts(j,:));
                if (r>0)
                    switch mod(kappa,2)
                        case 0
                            K(i,j) =   (r.^kappa)*(log(r));
                        case 1
                            K(i,j) = (r.^kappa);
                    end
                end
            end
        end
        
        [m,d] = size(landmarks);
        %U = zeros(m,n);
        U = lambda*ones(m,n);
        
        dUx = zeros(m,n);
        dUy = zeros(m,n);
        dUz = zeros(m,n);
        
        switch mod(kappa,2)
            case 0
                for i=1:m
                    for j=1:n
                        r = norm(landmarks(i,:) - ctrl_pts(j,:));
                        if (r>0)
                            U(i,j) =  (r.^kappa)*log(r);
                            dUx(i,j) = (r.^(kappa-2))*(landmarks(i,1)-ctrl_pts(j,1))*(kappa*log(r) +1);
                            dUy(i,j) = (r.^(kappa-2))*(landmarks(i,2)-ctrl_pts(j,2))*(kappa*log(r) +1);
                            if(d == 3)
                                dUz(i,j) = (r.^(kappa-2))*(landmarks(i,3)-ctrl_pts(j,3))*(kappa*log(r) +1);
                            end
                        end
                    end
                end
            case 1
                for i=1:m
                    for j=1:n
                        r = norm(landmarks(i,:) - ctrl_pts(j,:));
                        U(i,j) =   r.^kappa;
                        if(r > 0)
                            dUx(i,j) = (landmarks(i,1) - ctrl_pts(j,1))*r^(kappa-1)*kappa;
                            dUy(i,j) = (landmarks(i,2) - ctrl_pts(j,2))*r^(kappa-1)*kappa;
                            if(d == 3)
                                dUz(i,j) = (landmarks(i,3) - ctrl_pts(j,3))*r^(kappa-1)*kappa;
                            end
                        end
                    end
                end
        end
        
    case('clamped')
        %Work in progress
        switch d
            case 2
                for i=1:n
                    for j=1:n
                        a1 = norm(ctrl_pts(i,1:2))^2;
                        a2 = norm(ctrl_pts(j,1:2))^2;
                        a3 = ctrl_pts(i,1:2)*ctrl_pts(j,1:2)';
                        r = norm(ctrl_pts(i,1:2) - ctrl_pts(j,1:2));
                        if (r>0)
                            afull = sqrt(a1*a2 - 2 * a3+1)/r;
                            K(i,j) =  .5*(1 - r.^2 - r^2 * log(r^2));%r.^2 * (.5*(afull^2 - 1) - log(afull));
                        end
                    end
                end
            case 3
                for i=1:n
                    for j=1:n
                        a1 = norm(ctrl_pts(i,1:3))^2;
                        a2 = norm(ctrl_pts(j,1:3))^2;
                        a3 = ctrl_pts(i,1:3)*ctrl_pts(j,1:3)';
                        r = norm(ctrl_pts(i,1:3) - ctrl_pts(j,1:3));
                        if (r>0)
                            afull = sqrt(a1*a2 - 2 * a3+1)/r;
                            K(i,j) =  r.^2 * (.5*(afull^2 -1) + log(afull));
                        end
                    end
                end
        end
        
        if (nargin>=2)
            [m,d] = size(landmarks);
            %U = zeros(m,n);
            U = lambda*ones(m,n);
            
            dUx = zeros(m,n);
            dUy = zeros(m,n);
            dUz = zeros(m,n);
            
            switch d
                case 2
                    for i=1:m
                        thisX = landmarks(i,1:2);
                        for j=1:n
                            r = norm(thisX - ctrl_pts(j,1:2));
                            if (r>0)
                                a1 = norm(thisX)^2;
                                a2 = norm(ctrl_pts(j,1:2))^2;
                                a3 = landmarks(i,1:2)*ctrl_pts(j,1:2)';
                                afull = sqrt(a1*a2 - 2 * a3+1)/r;
                                U(i,j) =  .5*(1 - r.^2 - r^2 * log(r^2));%r.^2 * (.5*(afull^2 -1 ) - log(afull));
                                
                                Ax = ((thisX(1) * a2  - ctrl_pts(j,1)) * r / ...
                                    sqrt(a1*a2 - 2 * a3+1) - thisX(1)*afull) / ...
                                    (r^2 );
                                Ay = ( (2*thisX(2) * a2  -  ctrl_pts(j,2)) * r / ...
                                    sqrt(a1*a2 - 2 * a3+1) - thisX(2)*afull) / ...
                                    (r^2 );
                                dUx(i,j) = (thisX(1) - ctrl_pts(j,1)) * log(r^2);%(landmarks(i,1)-ctrl_pts(j,1))*(2*log(r) +1);
                                dUy(i,j) = (thisX(2) - ctrl_pts(j,2)) * log(r^2);%(landmarks(i,2)-ctrl_pts(j,2))*(2*log(r) +1);
                                
                                
                                %dUx(i,j) = 2 * (thisX(1) - ctrl_pts(j,1))* ...
                                %    (.5*(afull^2 -1) - log(afull)) + ...
                                %    r.^2 * (afull * Ax - Ax/afull);
                                
                                %dUy(i,j) = 2 * (thisX(2) - ctrl_pts(j,2))* ...
                                %    (.5*(afull^2 ) - log(afull)) + ...
                                %    r.^2 * (afull * Ay - Ay/afull);
                            end
                        end
                    end
                case 3
                    for i=1:m
                        for j=1:n
                            r = norm(landmarks(i,1:3) - ctrl_pts(j,1:3));
                            U(i,j) =   -r;
                            if(r > 0)
                                dUx(i,j) = - (landmarks(i,1) - ctrl_pts(j,1))/r;
                                dUy(i,j) = - (landmarks(i,2) - ctrl_pts(j,2))/r;
                                dUz(i,j) = - (landmarks(i,3) - ctrl_pts(j,3))/r;
                            end
                        end
                    end
            end
        end
    case('grbf')
        
        sigma = .35;
        denom=-2*sigma^2;
        [n, d]=size(landmarks); [m, d]=size(ctrl_pts);
        
        K0=repmat(ctrl_pts,[1 1 m])-permute(repmat(ctrl_pts,[1 1 m]),[3 2 1]);
        U0 = repmat(landmarks,[1 1 m])-permute(repmat(ctrl_pts,[1 1 n]),[3 2 1]);
        U=squeeze(sum(U0.^2,2))/denom;
        K = squeeze(sum(K0.^2,2))/denom;
        
        U=exp(U);
        K = exp(K);
        
        dUx = 2*squeeze(U0(:,1,:)).*U/denom;
        dUy = 2*squeeze(U0(:,2,:)).*U/denom;
        dUz =[];
        if(d == 3)
            dUz = 2*squeeze(U0(:,3,:)).*U/denom;
        end
        
end

