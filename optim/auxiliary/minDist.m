function d = minDist(X,lb,ub)
    X = (X-lb)./(ub-lb);
    PointDists = pdist(X);
    d = min(PointDists(:));
end