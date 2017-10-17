function gprMdl = fitGP(xTrain, yTrain, hyper, kernelFcn, optimizeNoiseHyper)
    doublings = 0;
    success = false;
    
    while ~success && doublings <= 20
        try
            noiseLB = hyper.noiseLB * 2 ^ doublings;
            if optimizeNoiseHyper
                gprMdl = fitrgp(xTrain, yTrain,...
                    'KernelParameters', [hyper.l, hyper.f],...
                    'Sigma', hyper.noise,...
                    'SigmaLowerBound', noiseLB,...
                    'KernelFunction', kernelFcn,...
                    'OptimizeHyperparameters', {'Sigma'},...
                    'HyperparameterOptimizationOptions', struct('ShowPlots',true,'Verbose',1, 'UseParallel', 1));
%                 global logNoise;
%                 logNoise = [logNoise; hyper.noise, gprMdl.ModelParameters.Sigma];
            else
                gprMdl = fitrgp(xTrain, yTrain,...
                    'KernelParameters', [hyper.l, hyper.f],...
                    'Sigma', hyper.noise,...
                    'SigmaLowerBound', noiseLB,...
                    'KernelFunction', kernelFcn); 
            end
    
        catch me
            doublings = doublings + 1;
        end
    end
    
    if ~success
        xxx = 1;
    end
end

