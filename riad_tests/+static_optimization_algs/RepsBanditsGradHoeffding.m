classdef RepsBanditsGradHoeffding
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['RepsBanditsGradHoeffding' '_' num2str(optimizerInput.confidence) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.initVar)];
        end
        
        function [perf, kls] = optimizeStruct(optimizerInput)
            [perf, kls] = static_optimization_algs.RepsBanditsGradHoeffding.optimize(optimizerInput.initDistrib, optimizerInput.confidence,...
                optimizerInput.nbSamplesPerIter, optimizerInput.maxIterReuse, optimizerInput.nbIter, optimizerInput.fun, optimizerInput.videoFile);
        end
        
        function [perf, kls] = optimize(initDistrib, confidence, nbSamplesPerIter, maxIterReuse, nbIter, fun, videoFile)
            perf = zeros(nbIter, 2);
            kls = zeros(nbIter, 1);
            nbSamplesForEval = 100;
            dim = initDistrib.getDim;
            lastPolicies{1} = initDistrib;
            allSamples = [];
            allVals = [];
            
            %% animation
            if(~isempty(videoFile))
                open(videoFile);
            end
            
            for iter = 1:nbIter
                %% evaluation of current distribution
                newSamples = lastPolicies{end}.getSamples(max(nbSamplesPerIter, nbSamplesForEval));
                vals = fun.eval(newSamples);
                newPerf = mean(vals);
                perf(iter, :) = [((iter-1) * nbSamplesPerIter) newPerf];
                
                %% new samples of the iteration
                newSamples = newSamples(1:nbSamplesPerIter, :);
                vals = vals(1:nbSamplesPerIter);
                allSamples = [allSamples; newSamples];
                allVals = [allVals; vals];
                %% current distrib plotting
                if(~isempty(videoFile))
                    frame = static_optimization_algs.RepsBanditsIProjection.getFrame(fun, lastPolicies{end}, allSamples);
                    writeVideo(videoFile, frame);
                end
                
                %% search distribution update                
                
                %%% importance sampling probas
                logInitProbas = lastPolicies{end}.getLogProbas(allSamples);
                logImportanceProbas = zeros(length(allVals), length(lastPolicies)); 
                for k = 1:length(lastPolicies)
                    logImportanceProbas(:, k) = lastPolicies{k}.getLogProbas(allSamples);
                end
                maxLogProbas = max(logImportanceProbas, [], 2);
                logImportanceProbas = log(sum(exp(bsxfun(@minus, logImportanceProbas, maxLogProbas)), 2));
                logImportanceProbas = bsxfun(@plus, logImportanceProbas, maxLogProbas);
                
                %%% initial configuration
                initPrecision = lastPolicies{end}.getPrecision;
                initMu = lastPolicies{end}.getMu;
                initDetPrecision = lastPolicies{end}.getDetPrecision;
                anderson0 = static_optimization_algs.ConcentrationInequalities.andersonLowerBound(allVals, logImportanceProbas, logInitProbas, confidence);
                                
                doNormalizeWeights = 1;
                
                %%% optimizing the I-Projection + Entropy + Confidence                
                optFun = @(eta)@(params) static_optimization_algs.IProjection.eval(params, initPrecision, initMu, allSamples, allVals, logImportanceProbas, eta, doNormalizeWeights);
                eta = 1e-5;
%                 else
%                     eta = eta / 2^20;
%                 end

                threshKlRet = 1e-3; 
                [muOptim, covOptim, detCovOptim, firstImprovement, eta] = static_optimization_algs.RepsBanditsGradHoeffding.loopOverEta(optFun, dim, initMu, initPrecision, initDetPrecision, anderson0, eta, confidence, allSamples, allVals, logImportanceProbas, threshKlRet);
                
                if(firstImprovement)
                    lastPolicies{end+1} = static_optimization_algs.Normal(muOptim, covOptim);
                    kls(iter) = static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                else
                    kls(iter) = 0;
                    lastPolicies{end+1} = lastPolicies{end};
                end
                fprintf('kl %g, eta %g, firstImprovement %d\n', kls(iter), eta, firstImprovement);

                %%% delete old policies
                if(length(lastPolicies) > maxIterReuse)
                    lastPolicies = lastPolicies(2:end);
                    allSamples = allSamples(nbSamplesPerIter+1:end, :);
                    allVals = allVals(nbSamplesPerIter+1:end, :);
                end
            end
            if(~isempty(videoFile))
                frame = static_optimization_algs.RepsBanditsIProjection.getFrame(fun, lastPolicies{end}, newSamples);
                writeVideo(videoFile, frame);
                close(videoFile);
            end
        end
        
        function [muOptimRet, covOptimRet, detCovOptimRet, firstImprovement, eta] = loopOverEta(fun, dim, initMu, initPrecision, initDetPrecision, anderson0, eta, confidence, newSamples, vals, logImportanceProbas, threshKlRet)
            optOptim =  optimset('GradObj','on', 'DerivativeCheck', 'off', 'Display', 'off');
            firstImprovement = false;
            muOptimRet = 0;
            covOptimRet = 0;
            detCovOptimRet = 0;
            params0 = [initMu'; initPrecision(tril(ones(dim)) == 1)];     
            klUp = inf; klLb = -inf;
            incumbent = anderson0;
            for iterEta = 1:40
                %% for given eta, optimize iprojection
                params = fminunc(fun(eta), params0, optOptim);
                muOptim = params(1:dim)';
                precisionOptim = zeros(dim);
                precisionOptim(tril(ones(dim)) == 1) = params(dim+1:end);
                precisionOptim = precisionOptim + precisionOptim' - diag(diag(precisionOptim));
                covOptim = inv(precisionOptim);
                detCovOptim = det(covOptim);
                
                %% compute lower bound of the sample based average reward estimate
                logCurrProbas = static_optimization_algs.Normal.getLogProbasPrecision(newSamples, muOptim, precisionOptim, 1/detCovOptim);
                lb = static_optimization_algs.ConcentrationInequalities.andersonLowerBound(vals, logImportanceProbas, logCurrProbas, confidence);
                
                %% saving best params if improvement
                if(lb > incumbent)
                    incumbent = lb;
                    etaUp = eta;
                    klLb = static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                    muOptimRet = muOptim;
                    covOptimRet = covOptim;
                    detCovOptimRet = detCovOptim;
                    if(~firstImprovement)
                        firstImprovement = true;
                        if(iterEta ~= 1)
                            etaLb = eta / 2;
                        else
                            etaLb = 0;
                        end
                    end
                elseif(firstImprovement)
                    etaLb = eta;
                    klUp = static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                end
                
                %% computing new eta
                if(~firstImprovement)
                    kl = static_optimization_algs.RepsBanditsGradReverse.KL(initPrecision, initDetPrecision, initMu, covOptim, detCovOptim, muOptim);
                    if(kl < threshKlRet)
                        return;
                    end
                    eta = eta * 2;
                else 
                    eta = (etaUp - etaLb) / 2 + etaLb;                
                end
                
                if(klUp < klLb)
                    warning('This should not happen');
                    return;
                elseif(klUp - klLb < threshKlRet)
                    return;
                end
            end
        end
    end
end
