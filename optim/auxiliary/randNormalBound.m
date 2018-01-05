function ret = randNormalBound(lb,ub,n)
    ret = rand(n,size(lb,2)) .* (ub - lb) + lb;
end