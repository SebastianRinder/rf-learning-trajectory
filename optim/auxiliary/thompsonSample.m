%% author: Sebastian Rinder

function val = thompsonSample(testX, knownX, knownY, trajectories, L, alpha, opts)
    persistent sampleX;
    persistent sampleY;
    
    knownX = [knownX; sampleX];
    knownY = [knownY; sampleY];
    
    [meanVec, covarianceVec] = gaussianProcess(testX, knownX, knownY, trajectories, L, alpha, opts);
    val = meanVec + sqrt(covarianceVec) .* randn;
        
    sampleX = [sampleX; testX];
    sampleY = [sampleY; val];
end