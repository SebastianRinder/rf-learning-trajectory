function [mean, variance] = gaussianProcessModel(testX, knownX, knownY, opts)
    Ds = opts.covarianceFcn(testX, knownX, opts)';
    Ks = scaleKernel(Ds, opts);
    
    if nargout == 1
        alpha = opts.alpha;
        mean = Ks' * alpha;      %Rasmussen, Williams 2005
        
%             K = opts.L;
%             Knoise = K + opts.hyperN.*eye(size(K));
%             mean = Ks' * (Knoise \ knownY);
    else
        Kss = scaleKernel(0, opts);
        L = opts.L;
        Lk = L \ Ks;       %Python example, Nando de Freitas 2013
        mean = Lk' * (L \ knownY);
        variance = (Kss - sum(Lk.^2))';

%             K = opts.L;
%             Knoise = K + opts.hyperN.*eye(size(K));
%             mean = Ks' * (Knoise \ knownY);
%             variance = Kss - diag(Ks' * (Knoise \ Ks));
    end
end