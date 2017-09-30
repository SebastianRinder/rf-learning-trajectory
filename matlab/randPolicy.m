function ret = randPolicy(lb,ub,n)
    dim = size(lb,2);
    ret = rand(n,dim) .* (ub - lb) + lb;
end