function ret = randTheta(lb,ub)
    dim = size(lb,2);
    ret = rand(1,dim) .* (ub - lb) + lb;
end