function [EI, meanMdl, stdMdl] = acquisition(X, trainX, trainY, minMeanMdl, covarianceFcn, K, hyper, opts)
    [meanMdl, stdMdl] = gaussianProcess(X, trainX, trainY, covarianceFcn, K, hyper, opts);
    
    stdv = sqrt(max(0, stdMdl.^2 - hyper.noise.^2));
    gamma = (minMeanMdl - meanMdl)./stdv;
    PI = normcdf(gamma, 0, 1);
    EI = stdv.*(gamma.*PI + normpdf(gamma, 0, 1));
end