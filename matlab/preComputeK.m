function [L, alpha] = preComputeK(knownX, knownY, opts)
%     noise = std(knownY)/sqrt(2);
    K = opts.covarianceFcn(knownX, knownX, opts);
    Knoise = K + opts.hyper(3).*eye(size(K));
    [L,p] = chol(Knoise, 'lower');
    if p > 0
        L = [];
    end
    
    if nargout > 1
        alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
    end
end