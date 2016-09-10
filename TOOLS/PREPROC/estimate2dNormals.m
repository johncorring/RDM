function [v T] = estimate2dNormals(x, k, ordered, smoothed, oriented)
if(~exist('smoothed', 'var'))
    smoothed = 0;
end
if(~exist('oriented', 'var'))
    oriented = 0;
end
[n d] = size(x);

x = x/max(max(abs(x)));

v = zeros(n,1);
T = zeros(n,d);

if(~ordered)
    %for each vector find its NNs and compute a normal vector (by ls's)
    for i = 1:n
        sqrds = sum((x - repmat(x(i,:), n, 1)).^2,2);
        [vals indxs] = sort(sqrds, 'ascend');
        
        k = min(k,size(x,1)-1);
        m = bsxfun(@times,sqrt(vals(2:k+1)), x(indxs(2:k+1),:)-repmat(x(i,:), k,1));
        [u s a] = svd(m'*m - eye(2));
        
        T(i,:) = u(1,:);
        v(i) = atan2(T(i,2),T(i,1));
    end
else
    if(~smoothed)
        %compute the normal to the subsequent vec
        for i = 2:n
            T(i,:) = x(i,:) - x(i-1,:);
            T(i,:) = T(i,:)/norm(T(i,:));
            v(i) = atan2(-T(i,1), T(i,2));
        end
        v(1) = (v(end)+v(2))/2;
    else
        %compute the normal to the subsequent vec
        x = [x(end:-1:end-1,:); x; x(1:2,:)];
        for i = 1:n
            j = i + 2;
            T(i,:) = mean( x([j-1 j j+1],:) - x([j-2 j- 1 j],:));
            T(i,:) = T(i,:)/norm(T(i,:));
            v(i) = atan2(-T(i,1), T(i,2));
        end
        v(1) = (v(end)+v(2))/2;
    end
end

if(oriented)
   v = reorient(x, k, v); 
end
end

function v = reorient(x, k, v)

figure;
A = squareform(pdist(x));
B = A-A;
for j = 1:size(x,1)
    
    [asorted inds] = sort(A(j,:),'ascend');
    B(j,:) = A(j,:) <= asorted(k);
    
end

newV = v;

Q = [1];
doneInds = [];

while( ~(isempty(Q)))
   p0 = x(Q(end),:);
   v0 = newV(Q(end));
   whereNeighbors = B(Q(end),:);
   doneInds = [doneInds; Q(end)];
   Q(end) = [];
   p0Inds = find(whereNeighbors);
   if(0)%showGraphics
       cla;   
       quiver(x(:,1), x(:,2), cos(newV), sin(newV), 'b');
       pause(.025)
   end
   for i = 1:length(p0Inds)
        if( ismember(p0Inds(i), doneInds) | ismember(p0Inds(i), Q))
           continue;
        end
        vout = alignNormals(newV(p0Inds(i)), v0, x(p0Inds(i),:), p0);
        newV(p0Inds(i)) = vout;
        Q = [Q; p0Inds(i)];
   end
end
v = newV

end

function vout = alignNormals(v1, v2, x1, x2);
    if(0)%showGraphics
        hold on;
        quiver([x1(1) x2(1)], [x1(2) x2(2)], cos([v1 v2]), sin([v1 v2]), 'r');
        pause(.025);
    end
    if(sum([cos(v1) sin(v1)] .* [cos(v2) sin(v2)]) > 0)%do for 3d
        vout = v1;
    else
        vout = v1 + pi;
    end
end