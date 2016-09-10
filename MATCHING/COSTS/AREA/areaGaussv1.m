function garea = areaGaussv1(A, B, scale)
%B is the moving template, A is the target, scale is the param

[n1 d] = size(A);
[n2 d] = size(B);
garea = 0;
for i = 1:n1
    for j = 1:n2
        
        %calc vecs
        vij = perp(B(j,:) - A(i,:));
        if(norm(vij) > 5*scale ) continue ; end
        
        mpij = (A(i,:) + B(j,:))/2;
        gij = exp( - sum((A(i,:) - B(j,:)).^2)/(4*scale^2));
        
        for k = 1:n1
            for l = 1:n2
                
                vkl = perp(B(l,:) - A(k,:));
                
                %calc midpoint
                mpkl = (A(k,:) + B(l,:))/2;
                mpijkl = (mpij + mpkl)/2;
                
                %calc coefficients
                gkl = exp( - sum((B(l,:) - A(k,:)).^2)/(4*scale^2));
                gijkl = exp( - sum((mpij - mpkl).^2)/(4*scale^2));
                
                %calc qf matrix
                Q = vkl' * vij;
                
                %calc integrals
                intQ = trace(Q)*scale^2 + mpijkl * (Q * mpijkl');
                intvij = vkl*B(l,:)' * (vij* (mpijkl)');
                intvkl = vij*B(j,:)' * (vkl* (mpijkl)');
                dijkl = vij*B(j,:)' * (vkl * B(l,:)');
                
                garea = garea + gij*gkl*gijkl*(intQ - intvij - intvkl + dijkl);
            end
        end
    end
end

garea = garea/(scale^8*n1*n2);


end

function pvec = perp(v)

pvec(1) = -v(2);
pvec(2) = v(1);

end