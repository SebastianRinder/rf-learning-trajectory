function [mean, variance] = gaussianProcessModel(testX, knownX, knownY, opts)
    Ds = opts.covarianceFcn(testX, knownX, opts)';
    Ks = scaleKernel(Ds, opts);
    
    if nargout == 1
        alpha = opts.alpha;
        try
            mean = Ks' * alpha;      %Rasmussen, Williams 2005
        catch me
            disp('gpm1');
        end
        
%             K = opts.L;
%             Knoise = K + opts.hyperN.*eye(size(K));
%             mean = Ks' * (Knoise \ knownY);
    else
        try
            Kss = scaleKernel(0, opts);
            L = opts.L;
            Lk = L \ Ks;       %Python example, Nando de Freitas 2013
            mean = Lk' * (L \ knownY);
            variance = (Kss - sum(Lk.^2))';
        catch me
            disp('gpm');
        end

%             K = opts.L;
%             Knoise = K + opts.hyperN.*eye(size(K));
%             mean = Ks' * (Knoise \ knownY);
%             variance = Kss - diag(Ks' * (Knoise \ Ks));
    end
end