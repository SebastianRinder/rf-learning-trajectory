classdef RepsBanditsGradNonParamIS
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsGradNonParamIS' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.RepsBanditsGradNonParamIS.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf, kls] = optimize(initDistrib, epsiKL, nbSamplesPerIter, nbInitSamples, nbIter, maxIterReuse, fun, videoFile)
            perf = zeros(nbIter, 2);
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 100;
            dim = initDistrib.getDim;
            currentDistrib = initDistrib;
            allSamples = [];
            allVals = [];
           
            %% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            for iter = 1:nbIter
                %% evaluation
                numberOfSamplesEval = max(nbSamplesPerIter, nbSamplesForEval);
                if(iter == 1)
                    numberOfSamplesEval = max(numberOfSamplesEval, nbInitSamples);
                end
                newSamples = currentDistrib.getSamples(numberOfSamplesEval);
                vals = fun.eval(newSamples);
                newPerf = mean(vals);
                
                %% samples
                if(iter == 1)
                    numberOfSamplesUpdate = nbInitSamples;
                    perf(iter, :) = [0 newPerf];
                else
                    numberOfSamplesUpdate = nbSamplesPerIter;
                    perf(iter, :) = [((iter-2) * nbSamplesPerIter + nbInitSamples) newPerf];
                end
                
                newSamples = newSamples(1:numberOfSamplesUpdate, :);
                vals = vals(1:numberOfSamplesUpdate);
                allSamples = [allSamples; newSamples];
                allVals = [allVals; vals];

                %% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.RepsBanditsIProjection.getFrame(fun, currentDistrib, allSamples);
                    writeVideo(videoFile, frame);
                end
                %% update
                if(~exist('eta', 'var'))
                    [eta, ~] = static_optimization_algs.RepsBandits.optimizeDual(allVals, epsiKL);
                end
                %%% optimizing the I-Projection
                initPrecision = currentDistrib.getPrecision;
                initDetPrecision = currentDistrib.getDetPrecision;
                initMu = currentDistrib.getMu;
                params0 = [initMu'; initPrecision(tril(ones(dim)) == 1)];
                
                %%%Importance weights
                logInitProbas = currentDistrib.getLogProbas(allSamples);
                logImportanceProbas = static_optimization_algs.NonParamImportanceSampling.getLogProbas(allSamples, initPrecision, logInitProbas);
%                 logImportanceProbas = ones(size(logInitProbas));

                optFun = @(eta)@(params) static_optimization_algs.IProjection.eval(params, initPrecision, initMu, allSamples, allVals, logImportanceProbas, eta);
                klFun = @(covOptim, detCovOptim, muOptim) static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                [~, eta, muOptim, covOptim, kl] = static_optimization_algs.RepsBanditsGrad.loopOverEta(optFun, klFun, dim, params0, eta, epsiKL, false);
                
                currentDistrib = static_optimization_algs.Normal(muOptim, covOptim);
                %%% delete old policies
                if(iter > maxIterReuse)
                    if(length(allVals) ~= maxIterReuse * nbSamplesPerIter)
                        allSamples = allSamples(nbInitSamples+1:end, :);
                        allVals = allVals(nbInitSamples+1:end, :);                        
                    else
                        allSamples = allSamples(nbSamplesPerIter+1:end, :);
                        allVals = allVals(nbSamplesPerIter+1:end, :);
                    end
                end
                kls(iter) = kl;
            end
            if(~isempty(videoFile))
                close(videoFile);
            end
        end
        
        
    end
end
