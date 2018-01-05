function [meanVec, covarianceVec, covarianceMat] = gaussianProcess(testX, knownX, trajectories, L, alpha, func)
    Ds = func.opts.distanceMat(testX, knownX, trajectories, false, func.opts)';
    Ks = func.opts.scaleKernel(Ds, func.opts.hyper);
    
    meanVec = Ks' * alpha; %Rasmussen, Williams 2005
    
    if nargout == 2
        Kss = func.opts.scaleKernel(0, func.opts.hyper);
        v = L \ Ks;       
        covarianceVec = (Kss - sum(v.^2))';    %Python example, Nando de Freitas 2013

    elseif nargout == 3
        covarianceVec = [];        
        Dss = func.opts.distanceMat(testX, testX, trajectories, true, opts)';
        Kss = func.opts.scaleKernel(Dss, func.opts.hyper);
        v = L \ Ks;
        covarianceMat = Kss - (v'*v); %Rasmussen, Williams 2005
        
    end
end