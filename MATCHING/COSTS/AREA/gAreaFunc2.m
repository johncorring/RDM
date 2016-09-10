function [areaGauss] = gAreaFunc2(A, B, scale)

[n1 d] = size(A);
[n2 d] = size(B);

%create 70x70 grid for matching
h = 50;
[xmax1] = max(A,[],1);
[xmax2] = max(B,[],1);
xmax = max(xmax1(1), xmax2(1));
ymax = max(xmax1(2), xmax2(2));

[xmin1] = min(A,[],1);
[xmin2] = min(B,[],1);
xmin = min(xmin1(1), xmin2(1));
ymin = min(xmin1(2), xmin2(2));

[x y] = meshgrid([xmin - (xmax - xmin)/4  : (xmax - xmin)/h  : xmax + (xmax + xmin)/4], ... 
    [ymin - (ymax - ymin)/4  : (ymax - ymin)/h  : ymax + (ymax + ymin)/4]);

areael = (xmax - xmin)/h   *  (ymax - ymin)/h  ;
areaGauss = x - x;

for i = 1:n1
    for j = 1:n2
        vix = -(y - A(i,2))/scale^2;
        viy = (x  - A(i,1))/scale^2;
        
        wjx = (x - B(j,1))/scale^2;
        wjy = (y - B(j,2))/scale^2;
        
        gFactorij = exp( - ((x - A(i,1)).^2 + (y - A(i,2)).^2)/(2*scale^2) );
        gFactorij = gFactorij .* exp( - ((x - B(j,1)).^2 + (y - B(j,2)).^2)/(2*scale^2) );
        
        areaGauss = areaGauss + (vix.*wjx + viy.*wjy).*gFactorij;
    end
end
areaGauss = sum(sum(areael*abs(areaGauss)/(n1*n2)));


end