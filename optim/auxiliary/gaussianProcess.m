function [meanVec, covarianceVec, covarianceMat] = gaussianProcess(testX, knownX, trajectories, L, alpha, func)
    Ds = func.opts.distanceMat(testX, knownX, trajectories, 'uk', func.opts)';
%     Ks = func.opts.scaleKernel(Ds, func.opts.hyper);
    Ks = func.opts.scaleKernel(Ds, [0,func.opts.hyperl.uk]);
    
    meanVec = Ks' * alpha; %Rasmussen, Williams 2005
    
    if nargout == 2
%         Kss = func.opts.scaleKernel(0, func.opts.hyper);
        Kss = func.opts.scaleKernel(0, [0,func.opts.hyperl.uk]);
        v = L \ Ks;       
        covarianceVec = (Kss - sum(v.^2))';    %Python example, Nando de Freitas 2013

    elseif nargout == 3
        covarianceVec = [];        
        Dss = func.opts.distanceMat(testX, testX, trajectories, 'uu', func.opts)';
        Kss = func.opts.scaleKernel(Dss, [0,func.opts.hyperl.uu]);
%         Kss = func.opts.scaleKernel(Dss, [0, 0]);
        v = L \ Ks;
        covarianceMat = Kss - (v'*v); %Rasmussen, Williams 2005
        
    end
end