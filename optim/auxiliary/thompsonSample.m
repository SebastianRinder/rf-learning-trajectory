%% author: Sebastian Rinder

function val = thompsonSample(X, knownX, knownY, trajectories, func)
    persistent sampleX;
    persistent sampleY;
        
    knownX = [knownX; sampleX];
    knownY = [knownY; sampleY];
    
    D = func.opts.distanceMat(knownX, knownX, trajectories, true, func.opts);
    K = func.opts.scaleKernel(D, func.opts.hyper);    
    [L, alpha] = getLowerCholesky(K, knownY, false, func.opts.sigmaNoiseSquared);
    
    [meanX, covarianceX] = gaussianProcess(X, knownX, trajectories, L, alpha, true, func);
    val = meanX + sqrt(covarianceX) .* randn;
        
    sampleX = [sampleX; X];
    sampleY = [sampleY; val];
end