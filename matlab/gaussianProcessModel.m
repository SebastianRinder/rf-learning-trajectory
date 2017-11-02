function [mean, variance] = gaussianProcessModel(testX, knownX, knownY, opts)
    Ks = opts.covarianceFcn(testX, knownX, opts)';
    
    if nargout == 1
        alpha = opts.alpha;
        mean = Ks' * alpha;      %Rasmussen, Williams 2005
    else
        L = opts.L;
        Kss = opts.covarianceFcn(testX(1,:), testX(1,:), opts);
        
        Lk = L \ Ks;       %Python example, Nando de Freitas 2013
        mean = Lk' * (L \ knownY);
        variance = (Kss - sum(Lk.^2))';
%         K = opts.covarianceFcn(knownX, knownX, opts);
%         mean = Ks * inv(K + 1e-6.*eye(size(K))) * knownY;
    end
end