function [meanVec, covarianceVec, covarianceMat] = gaussianProcess(testX, knownX, knownY, trajectories, L, alpha, opts)
    Ds = opts.distanceMat(testX, knownX, trajectories, false, opts)';
    Ks = opts.scaleKernel(Ds, opts.hyper);
    
    if nargout == 1
        meanVec = Ks' * alpha;      %Rasmussen, Williams 2005

    elseif nargout == 2
        Kss = opts.scaleKernel(0, opts.hyper);
        Lk = L \ Ks;       
        meanVec = Lk' * (L \ knownY);
        covarianceVec = (Kss - sum(Lk.^2))';    %Python example, Nando de Freitas 2013

    elseif nargout == 3
        meanVec = Ks' * alpha;        
        covarianceVec = [];        
        Dss = opts.distanceMat(testX, testX, trajectories, true, opts)';
        Kss = opts.scaleKernel(Dss, opts.hyper);
        v = L \ Ks;
        covarianceMat = Kss - (v'*v);     %GP_JPRED line 256
        
    end
end