function [newtemplate newVtemplate newVtarget energies] = solveDirectFreeFormL2(...
            template, Vtemplate, target, Vtarget, scale, lambda)

        %%iterate between computing the distance and solving for m, v
        %%start with uniform mixture ratios
        [n1 d] = size(template);
        [n2 d] = size(target);
        
        ptemp = ones(n1,1)/n1;
        ptarg = ones(n2,1)/n2;
        
        newptemp = ptemp;
        
        energies = [];
        energies0 = [];
        energies3 = [];
        energ1 = inf;
        distDiff = inf;
        
        newtemplate = template;
        newVtemplate = Vtemplate;
        
        newVtarget = Vtarget;
        
        newscale = scale;
        newlambda = lambda;
        
        templateMats = repmat(eye(2), n1, 1);
        targetMats = repmat(eye(2), n2, 1);
        
        for a = 1:n1
            
            templateMats(2*(a-1)+1,:) = [1 0]/newscale^2;
            templateMats(2*(a-1)+2,:) = [0 1]/newscale^2;
            
            %templateMats(2*(a-1)+1,:) = newVtemplate(a,:)/(norm(newVtemplate(a,:))*newscale^2)*2;
           % templateMats(2*(a-1)+2,:) = [-newVtemplate(a,2) newVtemplate(a,1)]/(norm(newVtemplate(a,:))*newscale^2);
        end
        for b = 1:n2
            targetMats(2*(b-1)+1,:) = [1 0]/newscale^2;
            targetMats(2*(b-1)+2,:) = [0 1]/newscale^2;
            
        end
        
        figure; 
        
        %TODO:tolerance param
        while(distDiff > 1e-6)
            %distance step
            damping = newscale/scale;
            [energ0 gradm0 gradv0 gradpb0] = gab_costfunct(newtemplate, target, newVtemplate,... 
                newVtarget,newptemp, ptarg, newscale , newlambda);
            
            [energ31 gradm31 gradv31 gradpb31] = gab_costfunct2(newtemplate, target, newVtemplate,... 
                newVtarget,templateMats, targetMats, newptemp, ptarg,  newlambda );
            
            [energ gradm gradv gradpb] = gab_costfunct3(newtemplate, target, newVtemplate,... 
                newVtarget,templateMats, targetMats, newptemp, ptarg,  newlambda );
            
            energies = [energies; energ];
            energies0 = [energies0; energ0];
            energies3 = [energies3; energ31];
            
            distDiff = abs(energ1 - energ)
            energ1 = energ;
            [energ20 gradm20 gradv20 gradpb20] =  gab_costfunct(target, newtemplate, newVtarget,... 
                newVtemplate, ptarg, newptemp, newscale ,newlambda);
            
            [energ32 gradm3 gradv3 gradpb3] =  gab_costfunct2(target, newtemplate, newVtarget,... 
                newVtemplate, targetMats, templateMats, ptarg, newptemp, newlambda );
            
            [energ2 gradm2 gradv2 gradpb2] =  gab_costfunct3(target, newtemplate, newVtarget,... 
                newVtemplate, targetMats, templateMats, ptarg, newptemp, newlambda );
            
            newVtarget = newVtarget - gradv2;
            newtemplate = newtemplate - gradm;
            newVtemplate = newVtemplate - gradv;

             
            
            %newptemp = newptemp - .0001*gradpb;
            %newptemp = newptemp/sum(newptemp);
            
            cla
            quiver(newtemplate(:,1), newtemplate(:,2), newVtemplate(:,1), ...
                newVtemplate(:,2), 'bx');
            hold on;
            quiver(target(:,1), target(:,2), newVtarget(:,1), newVtarget(:,2), ...
                'rs');
            
            for a = 1:n1
                [U S V] = svd(templateMats(2*(a-1)+1:2*a,:));
               hold on;
               quiver(newtemplate(a,1), newtemplate(a,2), U(1,1), U(1,2), .005*S(1,1));
               hold on;
               quiver(newtemplate(a,1), newtemplate(a,2), U(2,1), U(2,2), .005*S(2,2));
            end
            pause(.1);
            
            newscale = newscale * .99;
            newlambda = newlambda * .99;%sqrt(newscale);%newlambda * sqrt(.99);%
            
            
        for a = 1:n1
            
            templateMats(2*(a-1)+1,:) = [1 0]/newscale^2;
            templateMats(2*(a-1)+2,:) = [0 1]/newscale^2;
            
            %.1/.9 split can be a parameter
            Rot = [newVtemplate(a,:); [-newVtemplate(a,2) newVtemplate(a,1)]];
            templateMats(2*(a-1)+1:2*(a-1)+2,:) = Rot' * [.25/newscale^2 0 ; ...
                0 .75/newscale^2] * Rot;
            %templateMats(2*(a-1)+1,:) = newVtemplate(a,:)/(norm(newVtemplate(a,:))*newscale^2)*2;
            %templateMats(2*(a-1)+2,:) = [-newVtemplate(a,2) newVtemplate(a,1)]/(norm(newVtemplate(a,:))*newscale^2);
        end
        for b = 1:n2
            targetMats(2*(b-1)+1,:) = [1 0]/newscale^2;
            targetMats(2*(b-1)+2,:) = [0 1]/newscale^2;
            
        end
            newscale
            newlambda
        end%repeat
        
        
end

function f = ip2ptsGabor(ma, mb, va, vb, scale, lambda)

f = exp( - (norm(ma - mb)/(2*scale))^2 - ((scale/(2*lambda)) * norm(va - vb))^2 ... 
    + sqrt(-1) * sum((va + vb).*(ma - mb))/(2*lambda) );

end

function f = ip2ptsetsGabor(Ma, Mb, Va, Vb, scale, lambda)


        [n1 d] = size(Ma);
        [n2 d] = size(Mb);
        f = 0;
        
        for k = 1:n1
            for l = 1:n2
                f = f + ip2ptsGabor(Ma(k,:), Mb(l,:), Va(k,:), Vb(l,:), scale, lambda);
            end
        end

end


function [f g gv gpb] = gab_costfunct(A, B, Va, Vb, Pa, Pb, scale, lambda)
tic
    [f1 g1 gv1 gpb1] = GaborTransform2(A, A, Va, Va, Pa, Pa, lambda, scale);
    [f2 g2 gv2 gpb2] = GaborTransform2(A, B, Va, Vb, Pa, Pb, lambda, scale);
    [f3] = GaborTransform2(B, B, Vb, Vb, Pb, Pb, lambda, scale);
    toc
    f = f1 - 2*f2 + f3;
    g = 2*g1 - 2*g2;
    gv = 2*gv1 - 2*gv2;
    gpb = 2*gpb1 - 2*gpb2;
end

function [f g gv gpb] = gab_corr_costfunct(A, B, Va, Vb, Pa, Pb, scale, lambda)
    [f1 g1 gv1 gpb1] = GaborTransform2(A, A, Va, Va, Pa, Pa, lambda, scale);
    [f2 g2 gv2 gpb2] = GaborTransform2(A, B, Va, Vb, Pa, Pb, lambda, scale);
    [f3] = GaborTransform2(B, B, Vb, Vb, Pb, Pb, lambda, scale);
    f = -f2/((f1)*(f3));
    g = -(g2 *f1 * f3 - g1 * f2*f3)/(f1^2 * f3^2);%2*g1 - 2*g2;
    gv = -(gv2 *f1 * f3 - gv1 * f2*f3)/(f1^2 * f3^2);
    gpb = -(gpb2* f1 * f3 - gpb1 * f2*f3)/(f1^2 * f3^2);
end


function [f g gv gpb] = gab_costfunct2(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda)
tic
    [f1 g1 gv1 gpb1] = GaborTransform3(A, A, Va, Va, Ma, Ma, Pa, Pa, lambda);
    [f2 g2 gv2 gpb2] = GaborTransform3(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda);
    [f3 g3 gv3 gpb3] = GaborTransform3(B, B, Vb, Vb, Mb, Mb, Pb, Pb, lambda);
toc
    f = f1 - 2*f2 + f3;
    g = 2*g1 - 2*g2;
    gv = 2*gv1 - 2*gv2;
    gpb = 2*gpb1 - 2*gpb2;
end


function [f g gv gpb] = gab_costfunct3(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda)
tic
    [f1 g1 gv1 gpb1] = GaborTransform3_mat(A, A, Va, Va, Ma, Ma, Pa, Pa, lambda);
    [f2 g2 gv2 gpb2] = GaborTransform3_mat(A, B, Va, Vb, Ma, Mb, Pa, Pb, lambda);
    [f3 g3 gv3 gpb3] = GaborTransform3_mat(B, B, Vb, Vb, Mb, Mb, Pb, Pb, lambda);
    f1 
    f2 
    f3
toc
    f = f1 - 2*f2 + f3;
    g = 2*g1 - 2*g2;
    gv = 2*gv1 - 2*gv2;
    gpb = 2*gpb1 - 2*gpb2;
end
