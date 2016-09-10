function outpts = normalizePtSet(pts)
    [n d] = size(pts);
    %interpret as a flipped pointset
    if(d > n) 
        pts = pts'; 
    end
    
    outpts = pts - repmat(mean(pts,1), size(pts,1), 1);
    for k = 1:d
        outpts(:,k) = outpts(:,k)./max(abs(outpts(:,k)));
    end
end