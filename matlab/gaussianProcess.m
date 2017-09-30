function [meanPostY, stdPostY] = gaussianProcess(postX, priorX, priorY, covarianceFun, K, hyper, traj)
%     if isequal(covarianceFun, @sqExpCovariance)
%         Ks = covarianceFun(postX, priorX, [], hyper);
%         Kss = covarianceFun(postX, postX, [], hyper);
%     else
        Ks = covarianceFun(postX, priorX, hyper, traj);
        Kss = covarianceFun(postX, postX, hyper, traj);
%     end
    
    noise = eye(size(K)) .* hyper.noise;
    
    meanPostY = Ks * ((K + noise) \ priorY);
    varPostY = Kss - Ks * ((K + noise) \ Ks');
    stdPostY = sqrt(max(0,varPostY));
end