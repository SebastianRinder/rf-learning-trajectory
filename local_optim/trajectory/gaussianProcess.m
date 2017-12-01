function [meanVec, covarianceVec, covarianceMat] = gaussianProcess(testX, knownX, knownY, trajectories, opts)
    Ds = trajectoryCovariance(testX, knownX, trajectories, opts)';
    Ks = scaleKernel(Ds, opts);
    
    if nargout == 1
        alpha = opts.alpha;
        meanVec = Ks' * alpha;      %Rasmussen, Williams 2005

    elseif nargout == 2
        Kss = scaleKernel(0, opts);
        L = opts.L;
        Lk = L \ Ks;       %Python example, Nando de Freitas 2013
        meanVec = Lk' * (L \ knownY);
        covarianceVec = (Kss - sum(Lk.^2))';

    elseif nargout == 3
        alpha = opts.alpha;
        meanVec = Ks' * alpha;      %Rasmussen, Williams 2005
        
        covarianceVec = [];
        
        Dss = trajectoryCovariance(testX, testX, trajectories, opts)';
        Kss = scaleKernel(Dss, opts);
        v = opts.L \ Ks;
        covarianceMat = Kss - (v'*v);     %GP_JPRED line 256
        
    end
end