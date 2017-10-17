function [meanPostY, stdPostY] = gaussianProcess(postX, priorX, priorY, covarianceFcn, K, hyper, opts)
%     if isequal(covarianceFun, @sqExpCovariance)
%         Ks = covarianceFun(postX, priorX, [], hyper);
%         Kss = covarianceFun(postX, postX, [], hyper);
%     else
        Ks = covarianceFcn(postX, priorX, hyper, opts);
        Kss = covarianceFcn(postX, postX, hyper, opts);
%     end
    
    noise = eye(size(K)) .* hyper.noise;
    
    meanPostY = Ks * ((K + noise) \ priorY);
    if nargout > 1
        varPostY = Kss - Ks * ((K + noise) \ Ks');
        stdPostY = sqrt(max(0,varPostY));
    end
end