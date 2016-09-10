function [theta phi T] = est3dNormals(x, k)
    
    [n d] = size(x);
    
    %x = x - repmat(mean(x), n, 1);
    x = x/max(max(abs(x)));
    
    T = zeros(n,d);
    
    theta = zeros(n,1);
    phi = zeros(n,1);
    
    %for each vector find its NNs and compute a normal vector (by ls's)
    for i = 1:n
        sqrds = sum((x - repmat(x(i,:), n, 1)).^2,2);
        [vals indxs] = sort(sqrds, 'ascend');
        
        k = min(k,size(x,1)-1);
        m = bsxfun(@times,exp(-vals(2:k+1)), x(indxs(2:k+1),:)-repmat(x(i,:), k,1));
        m = m./sum(sum(m));
        %[u s a] = svd(m'*m - eye(d));
        [pc coeff scores] = princomp(m);%x(indxs(2:k+1),:)-repmat(x(i,:), k,1));
        %diag(s)
        T(i,:) = pc(:,end);%sum( repmat(pc(:,end), 1, k)' .* (x(indxs(2:k+1),:)-repmat(x(i,:), k,1)));%u(:,end);%
        [th phi2 r] = cart2sph(T(i,1), T(i,2), T(i,3));
        theta(i) = th;
        phi(i) = phi2;
    end


end