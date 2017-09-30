function gprMdl = fitGP(xTrain, yTrain, hyper, kernelFcn, optimizeNoiseHyper)
    if optimizeNoiseHyper
        gprMdl = fitrgp(xTrain, yTrain,...
            'KernelParameters', [hyper.l, hyper.f],...
            'Sigma', hyper.noise,...
            'SigmaLowerBound', hyper.noiseLB,...
            'KernelFunction', kernelFcn,...
            'OptimizeHyperparameters', {'Sigma'},...
            'HyperparameterOptimizationOptions', struct('ShowPlots',true,'Verbose',1, 'UseParallel', 1));
        global logNoise;
        logNoise = [logNoise; hyper.noise, gprMdl.ModelParameters.Sigma];
    else
        gprMdl = fitrgp(xTrain, yTrain,...
            'KernelParameters', [hyper.l, hyper.f],...
            'Sigma', hyper.noise,...
            'SigmaLowerBound', hyper.noiseLB,...
            'KernelFunction', kernelFcn); 
    end
end

