function [meanPostY, varPostY] = gaussianProcess(postX, priorY, priorX, covarianceFun, K, hyper, traj)
    if isequal(covarianceFun, @sqExpCovariance)
        Ks = covarianceFun(postX, priorX, [], hyper);
        %Kss = covarianceFun(postX, postX, hyper);
    else
        Ks = covarianceFun(postX, priorX, traj, hyper);
        %Kss = covarianceFun(postX, postX, traj, hyper);
    end
    
    Kss = 1;
    noise = eye(size(K)) .* 1e-6;
    
    meanPostY = Ks * (inv(K + noise) * priorY);
    varPostY = Kss - Ks * (inv(K + noise) * Ks');
    if varPostY<0
        varPostY = 0;
    end
end